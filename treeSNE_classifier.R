treeSNE_classifier <- function(runParams, classifierParams) {
	# treeSNE_classifier.R

	######################################################################################################
	# ADMIN
	######################################################################################################

	library(ggplot2)
	library(reshape2)
	library(data.table)
	library(RColorBrewer)
	library(stringr)
	library(penalized)
	library(nnet)
	library(glmnet)

	if (Sys.info()[['sysname']] == 'Darwin') {
		rootPath =  '~/Documents'; fileSep = '/'
	} else {
		rootPath = '\\\\d.ethz.ch\\dfs\\Groups\\biol\\sysbc\\users\\macnairw\\Documents'; fileSep = '\\'
	}

	set.seed(runParams$seed)


	######################################################################################################
	# FUNCTIONS
	######################################################################################################

	get_weights_outputs <- function(saveStem) {

		weightsFile = file.path(outputDir, paste0(saveStem, '_weights.txt'), fsep = fileSep);
		weightsOutputs = data.table(read.delim(weightsFile, header = TRUE, check.names = FALSE))
		return(weightsOutputs)
	}

	get_cluster_outputs <- function(saveStem) {
		clustersFile = file.path(outputDir, paste0(saveStem, '_clusters.txt'), fsep = fileSep);
		clustersOutputs = data.table(read.delim(clustersFile, header = TRUE, check.names = FALSE))
		return(clustersOutputs)

	}

	calc_coeffs_single_cluster <- function(clusterID, scaledClustersOutputs) {
		# make data binary
		responses = factor(scaledClustersOutputs$cluster == clusterID)
		penalizeds = as.matrix(scaledClustersOutputs[ ,-1, with = FALSE])
		
		cvfit = cv.glmnet(x = as.matrix(penalizeds), y = responses, nfolds = 10,
			family="binomial", alpha = 1, nlambda = 100, type.measure = classifierParams$type.measure) # type.multinomial=c("ungrouped","grouped")
		plot(cvfit)
		return(cvfit)
	}

	calc_fits <- function(clustersOutputs) {

		markers = clustersOutputs[ ,-1, with = FALSE]
		scaleFactors = list(
			min = apply(markers, 2, min),
			range = apply(markers, 2, function(x) quantile(x, probs = 0.95)) - apply(markers, 2, min)
			)
		scaledClustersOutputs = data.table(
			cluster = as.factor(clustersOutputs$cluster), 
			scale(markers, scaleFactors$min, scaleFactors$range)
			)

		clusters = sort(unique(scaledClustersOutputs$cluster))
		fits = lapply(clusters, function(x) {calc_coeffs_single_cluster(x, scaledClustersOutputs)} )
		return(fits)
	}

	extract_coeffs <- function(fits) {
		coeffs = sapply(fits, function(x) { as.matrix(coef(x)) })
		rownames(coeffs) = c('Intercept', names(clustersOutputs)[-1])
		colnames(coeffs) = 1:ncol(coeffs)
		return(coeffs)
	}

	extract_errors <- function(fits) {
		errors = data.frame(misclass = sapply(fits, function(x) { min(x$cvm) }), cluster = 1:length(fits))
		return(errors)
	}

	prep_meltDT <- function(coeffs) {
		names(dimnames(coeffs)) = c('marker', 'clusterID')
		meltDT = data.table(melt(coeffs))

		# cluster non-zeros together
		binaryDist = dist(t(coeffs), method = 'binary')
		coeffClust = hclust(binaryDist)
		# change levels accordingly
		colnames(coeffs) = 1:ncol(coeffs)
		newClusterLabels = as.character(coeffClust$order)
		meltDT$clusterID = factor(meltDT$clusterID, levels = newClusterLabels)

		nonZeroCount = rowSums(coeffs!=0)
		meltDT$marker = factor(meltDT$marker, levels = meltDT$marker[order(nonZeroCount)])
		
		return(meltDT)
	}

	prep_meltDT_unsorted <- function(coeffs) {

		names(dimnames(coeffs)) = c('marker', 'clusterID')
		meltDT = data.table(melt(coeffs))

		nonZeroCount = rowSums(coeffs!=0)
		meltDT$marker = factor(meltDT$marker, levels = meltDT$marker[order(nonZeroCount)])
		
		return(meltDT)
	}

	plot_coeffs <- function(coeffs) {

		# prepare data
		meltDT = prep_meltDT(coeffs[-1,])

		# plot parameters
		# fillValues = brewer.pal(9, 'Spectral')
		fillValues = brewer.pal(11, 'RdGy')
		limitValue = max(abs(meltDT$value))
		plotH = 8; plotW = plotH;

		# plot
		if (dev.cur()!=1) {
			dev.off()
		}
			
		saveFilename = file.path(outputDir, paste0(saveStem, '_classifier', classifierParams$fileType), fsep = fileSep)
		g = ggplot(meltDT, aes(y = marker, x = as.factor(clusterID), fill = value)) +
			geom_tile() +
			scale_fill_gradientn(colours = fillValues, limits = c(-limitValue, limitValue)) +
			theme_minimal() +
			labs(list(y = NULL, x = 'clusterID', fill = 'L1-penalized\ncoefficient'))
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_coeffs_unsorted <- function(coeffs) {

		# prepare data
		meltDT = prep_meltDT_unsorted(coeffs[-1,])
		meltDT$marker = factor(meltDT$marker, levels = sort(levels(meltDT$marker)))

		# plot parameters
		# fillValues = brewer.pal(9, 'Spectral')
		fillValues = brewer.pal(11, 'RdGy')
		limitValue = max(abs(meltDT$value))
		plotH = 8; plotW = plotH;

		# plot
		if (dev.cur()!=1) {
			dev.off()
		}
		saveFilename = file.path(outputDir, paste0(saveStem, '_classifier_unsorted', classifierParams$fileType), fsep = fileSep)
		g = ggplot(meltDT, aes(y = marker, x = as.factor(clusterID), fill = value)) +
			geom_tile() +
			scale_fill_gradientn(colours = fillValues, limits = c(-limitValue, limitValue)) +
			theme_minimal() +
			labs(list(y = NULL, x = 'clusterID', fill = 'L1-penalized\ncoefficient'))
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_markers <- function(coeffs) {

		noIntercept = coeffs[-1,]
		# prepare data
		counts = data.frame(marker = rownames(noIntercept), nonZeros = rowSums(noIntercept!=0))
		counts$marker = factor(counts$marker, levels = counts$marker[order(counts$nonZeros, decreasing = TRUE)])

		# plot parameters
		plotH = 4; plotW = plotH*2;
		# plot
		if (dev.cur()!=1) {
			dev.off()
		}
		saveFilename = file.path(outputDir, paste0(saveStem, '_useful_markers', classifierParams$fileType), fsep = fileSep)
		g = ggplot(counts, aes(y = nonZeros, x = marker)) +
			geom_bar(stat='identity') +
			theme(
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
				) +
			labs(list(x = 'Marker', y = 'Total observed non-zero coefficients'))
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_errors <- function(errors) {

		# prep data
		naiveError = 1/nrow(errors)
		maxY = ceiling(100*max(naiveError, errors$misclass))/100

		# plot parameters
		plotH = 4; plotW = plotH*1.7;

		# plot
		if (dev.cur()!=1) {
			dev.off()
		}
		saveFilename = file.path(outputDir, paste0(saveStem, '_classifier_errors', classifierParams$fileType), fsep = fileSep)
		g = ggplot(errors, aes(y = misclass, x = as.factor(cluster))) +
			geom_bar(stat = 'identity') +
			geom_hline(yintercept = naiveError) +
			annotate("text", label = "Naive classifier error rate", x = quantile(errors$cluster, 0.25), y = naiveError*1.05, size = 3) +
			coord_cartesian(ylim = c(0,maxY)) +
			labs(list(y = 'Misclassification error', x = 'clusterID'))
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	save_coeffs <- function(coeffs) {
		outputFile = file.path(outputDir, paste0(saveStem, '_coeffs.txt'), fsep = fileSep)
		write.table(coeffs, file = outputFile, append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
	}


	######################################################################################################
	# PARAMETERS
	######################################################################################################
	
	outputDir = runParams$outputDir
	saveStem = runParams$saveStem


	######################################################################################################
	# CODE
	######################################################################################################
	
	# import data
	clustersOutputs = get_cluster_outputs(saveStem)
	# clustersOutputs = clustersOutputs[cluster %in% 1:4, 1:5, with = FALSE]
	weightsOutputs = get_weights_outputs(saveStem)

	# calculate regularized coefficients
	fits = calc_fits(clustersOutputs)
	coeffs = extract_coeffs(fits)
	errors = extract_errors(fits)

	# plots
	save_coeffs(coeffs)
	plot_coeffs(coeffs)
	plot_coeffs_unsorted(coeffs)
	plot_markers(coeffs)
	plot_errors(errors)

}


# runParams = list(
# 	saveStem = 'bodenmiller_test_run_phenotypic_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_test_run_phenotypic'
# 	)
# classifierParams = list(
# 	seed = 0,
# 	fileType = '.png', 
# 	nfolds = 10, 
# 	alpha = 1, 
# 	nlambda = 100, 
# 	type.multinomial="ungrouped", 
# 	type.measure="class"
# 	)
# treeSNE_classifier(runParams, classifierParams)

runParams = list(
	saveStem = 'prostate_run_split_large_mean',
	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller_PC\\prostate_run_split_large\\partial outputs'
	)
classifierParams = list(
	seed = 0,
	fileType = '.png', 
	nfolds = 10, 
	alpha = 1, 
	nlambda = 100, 
	type.multinomial="ungrouped", 
	type.measure="class"
	)
treeSNE_classifier(runParams, classifierParams)
