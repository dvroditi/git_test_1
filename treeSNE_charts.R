treeSNE_charts <- function(inputParameters) {
	# treeSNE_charts.R

	######################################################################################################
	# ADMIN
	######################################################################################################

	require(ggplot2)
	require(reshape2)
	require(data.table)
	require(RColorBrewer)
	require(stringr)

	if (Sys.info()[['sysname']] == 'Darwin') {
		rootPath =  '~/Documents'; fileSep = '/'
	} else {
		rootPath = '\\\\d.ethz.ch\\dfs\\Groups\\biol\\sysbc\\users\\macnairw\\Documents'; fileSep = '\\'
	}


	######################################################################################################
	# FUNCTIONS
	######################################################################################################

	get_tsne_outputs <- function(saveStem) {

		tsneFile = file.path(outputDir, paste0(saveStem, '_tsne.txt'), fsep = fileSep);
		tsneOutputs = data.table(read.delim(tsneFile, header = TRUE))
		return(tsneOutputs)
	}

	get_marker_outputs <- function(saveStem, markerType) {

		markerFile = file.path(outputDir, paste0(saveStem, '_', markerType, '_markers.txt'), fsep = fileSep);
		scaledMarkers = tryCatch(
			{
				markerOutputs = data.table(read.delim(markerFile, header = TRUE))
				scaledMarkers = data.table(scale(markerOutputs, apply(markerOutputs, 2, min), apply(markerOutputs, 2, function(x) quantile(x, probs = 0.95)) - apply(markerOutputs, 2, min)))
			},
			error=function(cond) {
				scaledMarkers = data.table()
				return(scaledMarkers)
			}
		)
		return(scaledMarkers)
	}

	get_weights_outputs <- function(saveStem) {

		weightsFile = file.path(outputDir, paste0(saveStem, '_weights.txt'), fsep = fileSep);
		weightsOutputs = data.table(read.delim(weightsFile, header = TRUE, check.names = FALSE))
		return(weightsOutputs)
	}

	calc_max_weights <- function(weightsOutputs) {
		sampleNames = colnames(weightsOutputs)
		max_idx <- function(row) {
			max_idx = sample(which(row == max(row)), 1)
			return(max_idx)
		}
		maxSample = apply(weightsOutputs, 1, function(row) {sampleNames[max_idx(row)]})
		maxValue = apply(weightsOutputs, 1, function(row) {max(row)[1]})
		maxWeights = data.table(weight = maxValue, sample = maxSample)
		return(maxWeights)
	}

	calc_hist <- function(weightsOutputs) {
		mode = apply(weightsOutputs, 1, function(row) { max(row) / sum(row) })
		noSamples = apply(weightsOutputs, 1, function(row) { sum(row > 0) })
		DT = data.table(mode, noSamples)
		return(DT)
	}

	plot_simple <- function(tsneOutputs, weightsOutputs) {
		
		totalWeights = rowSums(weightsOutputs)

		DT = cbind(tsneOutputs, Cells = totalWeights)
		# set up plotting variables
		# plotH = 8; plotW = plotH*1.6
		plotH = 8; plotW = plotH*0.8

		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_simple', '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = Y1, y = Y2, size = Cells, alpha = Cells)) +
			geom_point() +
			scale_size(range = c(0,12)) +
			xlab(NULL) + ylab(NULL) +
			# labs(list(colour = "Sample name", shape = "Sample category")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey() +
			theme(legend.position = "bottom")
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_by_sample <- function(tsneOutputs, weightsOutputs) {
		
		# set up data
		DT = cbind(tsneOutputs, weightsOutputs)
		meltDT = data.table(melt(DT, id = c('Y1', 'Y2'), variable.name = 'sample', value.name = 'weight'))
		# meltDT = meltDT[weight != 0]
		meltDT$sample = str_replace_all(meltDT$sample, '\\.', '-')
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		colourValues = rev(brewer.pal(9, 'Spectral'))
		
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_by_sample', '.png'), fsep = fileSep)
		g = ggplot(meltDT, aes(x = Y1, y = Y2, size = weight, colour = weight)) +
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
			geom_point() +
			facet_wrap( ~ sample, ncol = 5) +
			scale_size(range = c(0,10)) +
			scale_colour_gradientn(colours = colourValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Approx\nno. of\ncells", size = "Approx\nno. of\ncells", alpha = NULL)) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_by_sample2 <- function(tsneOutputs, weightsOutputs) {
		
		# set up data
		DT = cbind(tsneOutputs, weightsOutputs)
		meltDT = data.table(melt(DT, id = c('Y1', 'Y2'), variable.name = 'sample', value.name = 'weight'))
		# meltDT = meltDT[weight != 0]
		meltDT$sample = str_replace_all(meltDT$sample, '\\.', '-')
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		colourValues = rev(brewer.pal(9, 'Spectral'))
		
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_by_sample2', '.eps'), fsep = fileSep)
		g = ggplot(meltDT, aes(x = Y1, y = Y2, colour = weight)) +
			geom_point() +
			# facet_wrap( ~ sample, ncol = 5) +
			# scale_size(range = c(0,10)) +
			scale_colour_gradientn(colours = colourValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Proportional\ncell density")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_by_sample_together <- function(tsneOutputs, weightsOutputs) {
		
		# set up data
		maxWeights = calc_max_weights(weightsOutputs)
		DT = cbind(tsneOutputs, maxWeights)
		# DT$sample = str_replace_all(DT$sample, '\\.', '-')
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		nSamples = length(unique(DT$sample))
		colourValues = colorRampPalette(brewer.pal(8, 'Dark2'))(nSamples)
		shapeValues = rep(0:3, times = nSamples)[1:nSamples]
		
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_by_sample_together', '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = Y1, y = Y2, size = weight, colour = sample, shape = sample)) +
			geom_point() +
			scale_size(range = c(0,10)) +
			scale_colour_manual(values = colourValues) +
			scale_shape_manual(values = shapeValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Sample\nname", shape = "Sample\nname", size = "Approx\nno. of\ncells")) +
			ggtitle('Most common sample at each point') +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_by_marker <- function(tsneOutputs, weightsOutputs, markerOutputs, markerType) {
		out = tryCatch( 
			{
				# set up data
				totalWeights = rowSums(weightsOutputs)
				DT = cbind(tsneOutputs, markerOutputs, weight = totalWeights)
				meltDT = data.table(melt(DT, id = c('Y1', 'Y2', 'weight'), variable.name = 'marker', value.name = 'value'))
				# meltDT = meltDT[weight != 0]
				# meltDT$marker = str_replace_all(meltDT$marker, '\\.', '-')
				
				# set up plotting variables
				nMarkers = dim(markerOutputs)[2]
				nCol = ceiling(sqrt(nMarkers*1.5))
				plotW = 2*nCol; plotH = plotW/1.5
				colourValues = rev(brewer.pal(9, 'Spectral'))
				
				# plot
				saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_by_', markerType, '_marker', '.png'), fsep = fileSep)
				# g = ggplot(meltDT, aes(x = Y1, y = Y2, size = weight, colour = value)) +
				g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, colour = value, size = weight)) +
				# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
					geom_point() +
					facet_wrap( ~ marker, ncol = nCol) +
					scale_size(range = c(0,4)) +
					# scale_alpha_manual(guide = 'none') +
					scale_colour_gradientn(colours = colourValues, limits = c(0,1)) +
					xlab(NULL) + ylab(NULL) +
					labs(list(colour = "Marker\nvalue", size = "Approx\nno. of\ncells", alpha = "Approx\nno. of\ncells")) +
					# ggtitle(paste0(fileStem, ' transform ', ii)) +
					theme_grey()
				ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)
			},
			error=function(cond) {
				message(paste0('Did not save extra markers chart for ', saveStem, '; extra markers file empty.'))
				return(NA)
			}
		)

	}

	plot_no_samples <- function(tsneOutputs, weightsOutputs) {
		
		# set up data
		DT = cbind(tsneOutputs, calc_hist(weightsOutputs))
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		nCols = length(unique(DT$noSamples))
		colourValues = rev(colorRampPalette(brewer.pal(11, 'Spectral'))(nCols))
		
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_no_samples', '.png'), fsep = fileSep)
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, size = weight, colour = value)) +
		g = ggplot(DT, aes(x = Y1, y = Y2, colour = as.factor(noSamples), size = noSamples)) +
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
			geom_point(size = 3) +
			scale_colour_manual(values = colourValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "No. of\nsamples")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		colourValues = rev(brewer.pal(9, 'Spectral'))
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_tsne_mode', '.png'), fsep = fileSep)
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, size = weight, colour = value)) +
		g = ggplot(DT, aes(x = Y1, y = Y2, colour = mode)) +
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
			geom_point(size = 3) +
			scale_colour_gradientn(colours = rev(colourValues)) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Highest\nsample\nproportion")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	plot_histograms <- function(weightsOutputs) {
		
		# set up data
		DT = calc_hist(weightsOutputs)
		
		# set up plotting variables
		plotH = 5; plotW = plotH*1.6
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_hist_no_samples', '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = noSamples)) +
			geom_histogram(aes(y = ..count.. / sum(..count..)), binwidth = 1) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			xlab('Number of samples observed at cluster') + ylab('Proportion') +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

		# set up plotting variables
		plotH = 5; plotW = plotH*1.6
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_hist_mode', '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = mode)) +
			geom_histogram(aes(y = ..count.. / sum(..count..)), binwidth = 0.1) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			xlab('Highest proportion of cluster from one sample') + ylab('Proportion') +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}



	######################################################################################################
	# CODE
	######################################################################################################
	
	saveStem = inputParameters$saveStem
	outputDir = inputParameters$outputDir

	# import data
	tsneOutputs = get_tsne_outputs(saveStem)
	scaledUsedMarkers = get_marker_outputs(saveStem, 'used')
	scaledExtraMarkers = get_marker_outputs(saveStem, 'extra')
	weightsOutputs = get_weights_outputs(saveStem)

	# plots
	plot_simple(tsneOutputs, weightsOutputs)
	plot_by_sample(tsneOutputs, weightsOutputs)
	plot_by_sample2(tsneOutputs, weightsOutputs)
	plot_by_marker(tsneOutputs, weightsOutputs, scaledUsedMarkers, 'used')
	plot_by_marker(tsneOutputs, weightsOutputs, scaledExtraMarkers, 'extra')
	plot_by_sample_together(tsneOutputs, weightsOutputs)
	plot_no_samples(tsneOutputs, weightsOutputs)
	plot_histograms(weightsOutputs)

}


# inputParameters = list(
# 	saveStem = 'treeSNE_robust_distance_phenotype_only_mean',
# 	outputDir = paste0('\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Data\\Cytobank\\', 'treeSNE_robust_distance_phenotype_only'),
# 	experimentNo = 2
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'treeSNE_robust_distance_phenotype_only_1234_mean',
# 	outputDir = paste0('\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Data\\Cytobank\\', 'treeSNE_robust_distance_phenotype_only_1234'),
# 	experimentNo = 2
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = '10dim_1e4runs_treeSNE_mean',
# 	outputDir = '\\\\d.ethz.ch\\dfs\\Groups\\biol\\sysbc\\users\\macnairw\\Documents\\Outputs\\02_Cell_differentiation_simulation\\treeSNE_positive_control\\10dim_1e4runs_treeSNE'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'bodenmiller_test_run_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_test_run'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'bodenmiller_test_run_phenotypic_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_test_run_phenotypic'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = '10dim_1e4runs_treeSNE_dijkstra_test_mean',
# 	outputDir = '\\\\d.ethz.ch\\dfs\\Groups\\biol\\sysbc\\users\\macnairw\\Documents\\Outputs\\02_Cell_differentiation_simulation\\treeSNE_positive_control\\10dim_1e4runs_treeSNE'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'bodenmiller_test_run_functional_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_test_run_functional'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'bodenmiller_cd4_sunitinib_IFNa_concs_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_cd4_sunitinib_IFNa_concs'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'treeSNE_test_prostate_2_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmillier_PC\\treeSNE_test_prostate_2'
# 	)
# treeSNE_charts(inputParameters)

# inputParameters = list(
# 	saveStem = 'treeSNE_test_prostate_3_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmillier_PC\\treeSNE_test_prostate_3'
# 	)
# treeSNE_charts(inputParameters)

inputParameters = list(
	saveStem = 'treeSNE_test_prostate_3_mean',
	outputDir = '/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/treeSNE_test_prostate_3'
	)
treeSNE_charts(inputParameters)

