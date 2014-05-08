treeSNE_clusterer <- function(runParams, clusterParams) {
	# treeSNE_clusterer.R

	######################################################################################################
	# ADMIN
	######################################################################################################

	library(ggplot2)
	library(reshape2)
	library(data.table)
	library(RColorBrewer)
	library(stringr)
	library(fpc)

	if (Sys.info()[['sysname']] == 'Darwin') {
		rootPath =  '~/Documents'; fileSep = '/'
	} else {
		rootPath = '\\\\d.ethz.ch\\dfs\\Groups\\biol\\sysbc\\users\\macnairw\\Documents'; fileSep = '\\'
	}

	set.seed(clusterParams$seed)


	######################################################################################################
	# FUNCTIONS
	######################################################################################################

	get_tsne_outputs <- function(saveStem) {

		tsneFile = file.path(outputDir, paste0(saveStem, '_tsne.txt'), fsep = fileSep);
		tsneOutputs = data.table(read.delim(tsneFile, header = TRUE))
		return(tsneOutputs)
	}

	# get_marker_outputs <- function(saveStem, markerType) {

	# 	markerFile = file.path(outputDir, paste0(saveStem, '_', markerType, '_markers.txt'), fsep = fileSep);
	# 	markerOutputs = data.table(read.delim(markerFile, header = TRUE))
	# 	scaledOutputs = data.table(scale(markerOutputs, apply(markerOutputs, 2, min), apply(markerOutputs, 2, function(x) quantile(x, probs = 0.95)) - apply(markerOutputs, 2, min)))


	# 	return(scaledOutputs)
	# }

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

	get_combined_marker_outputs <- function(saveStem) {
		usedMarkers = get_marker_outputs(saveStem, 'used')
		extraMarkers = get_marker_outputs(saveStem, 'extra')
		if (dim(extraMarkers)[1] == 0) {
			markerOutputs = usedMarkers
		} else {
			markerOutputs = cbind(usedMarkers, extraMarkers)
		}
		return(markerOutputs)
	}

	get_weights_outputs <- function(saveStem) {

		weightsFile = file.path(outputDir, paste0(saveStem, '_weights.txt'), fsep = fileSep);
		weightsOutputs = data.table(read.delim(weightsFile, header = TRUE, check.names = FALSE))
		return(weightsOutputs)
	}

	calc_clusters <- function(tsneOutputs) {
		clusterApproach = clusterParams$approach
		switch (clusterApproach,
			'kmeans' = {
				k = clusterParams$k
				nstarts = clusterParams$nstarts
				clusters = kmeans(tsneOutputs, k, iter.max = 10, nstart = nstarts)
				return(clusters[['cluster']])
			},
			'pam' = {
				krange = clusterParams$krange
				usepam = (dim(tsneOutputs)[1] <= 2000)
				clusters = pamk(tsneOutputs, krange=krange, criterion="asw", usepam=usepam, # usepam best for small (<2000) datasets
							scaling=FALSE, critout = TRUE)
				return(clusters$pamobject$clustering)
			})
	}

	calc_clusters_k <- function(tsneOutputs, k) {
		clusterApproach = clusterParams$approach
		switch (clusterApproach,
			'kmeans' = {
				k = k
				nstarts = clusterParams$nstarts
				clusters = kmeans(tsneOutputs, k, iter.max = 10, nstart = nstarts)
				return(clusters[['cluster']])
			},
			'pam' = {
				krange = clusterParams$krange
				usepam = (dim(tsneOutputs)[1] <= 2000)
				clusters = pamk(tsneOutputs, krange=krange, criterion="asw", usepam=usepam, # usepam best for small (<2000) datasets
							scaling=FALSE, critout = TRUE)
				return(clusters$pamobject$clustering)
			})
	}

	plot_clusters <- function(tsneOutputs, clusters) {
		
		# set up data
		DT = cbind(tsneOutputs, cluster = clusters)
		nClusters = length(unique(clusters))
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		colourValues = colorRampPalette(brewer.pal(8, 'Dark2'))(nClusters)
		shapeValues = rep(0:3, times = nClusters)[1:nClusters]

		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_clusters', '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = Y1, y = Y2, colour = as.factor(cluster), shape = as.factor(cluster))) +
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
			geom_point() +
			scale_colour_manual(values = colourValues) +
			scale_shape_manual(values = shapeValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Cluster\nID", shape = "Cluster\nID")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)			

	}

	plot_clusters_k <- function(tsneOutputs, clusters, kk) {
		
		# set up data
		DT = cbind(tsneOutputs, cluster = clusters)
		nClusters = length(unique(clusters))
		
		# set up plotting variables
		plotH = 8; plotW = plotH*1.6
		colourValues = colorRampPalette(brewer.pal(8, 'Dark2'))(nClusters)
		shapeValues = rep(0:3, times = nClusters)[1:nClusters]

		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_clusters_', kk, '.png'), fsep = fileSep)
		g = ggplot(DT, aes(x = Y1, y = Y2, colour = as.factor(cluster), shape = as.factor(cluster))) +
		# g = ggplot(meltDT, aes(x = Y1, y = Y2, alpha = weight, size = weight, colour = weight)) +
			geom_point() +
			scale_colour_manual(values = colourValues) +
			scale_shape_manual(values = shapeValues) +
			xlab(NULL) + ylab(NULL) +
			labs(list(colour = "Cluster\nID", shape = "Cluster\nID")) +
			# ggtitle(paste0(fileStem, ' transform ', ii)) +
			theme_grey()
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)			

	}

	plot_clusters_proportions <- function(clustersOutputs, weightsOutputs) {

		# set up data
		DT = data.table(cluster = clustersOutputs, weightsOutputs)
		meltDT = data.table(melt(DT, id = 'cluster', variable.name = 'sample'))

		# plot parameters
		plotH = 8; plotW = plotH*1.3;
		nSamples = length(unique(meltDT$sample))
		colourValues = colorRampPalette(brewer.pal(8, 'Dark2'))(nSamples)
		# plot
		saveFilename = file.path(outputDir, paste0(saveStem, '_cluster_proportions.png'), fsep = fileSep)
		g = ggplot(meltDT, aes(y = value, x = as.factor(cluster), fill = factor(sample))) +
			geom_bar(stat = 'identity', position = 'fill') +
			scale_fill_manual(values = colourValues) +
			theme(
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
				) +
			labs(list(x = 'clusterID', y = 'Proportion of cluster', fill = 'Sample'))
		ggsave(saveFilename, g, height = plotH, width = plotW, dpi = 500)

	}

	save_clusters <- function(clusters, markerOutputs) {
		output = data.frame(cluster = clusters, markerOutputs)
		outputFile = file.path(outputDir, paste0(saveStem, '_clusters.txt'), fsep = fileSep)
		write.table(output, file = outputFile, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

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
	tsneOutputs = get_tsne_outputs(saveStem)
	markerOutputs = get_combined_marker_outputs(saveStem)
	weightsOutputs = get_weights_outputs(saveStem)

	if (clusterParams$approach == 'kmeans' & "krange" %in% names(clusterParams)) {
		for (k in clusterParams$krange) {
			# find clusters
			clustersOutputs = calc_clusters_k(tsneOutputs, k)
			# plots
			plot_clusters_k(tsneOutputs, clustersOutputs, k)
		}
	}
	else {
		# find clusters
		clustersOutputs = calc_clusters(tsneOutputs)

		# plots
		plot_clusters(tsneOutputs, clustersOutputs)
		plot_clusters_proportions(clustersOutputs, weightsOutputs)
		save_clusters(clustersOutputs, markerOutputs)
	}
}



# runParams = list(
# 	saveStem = 'bodenmiller_test_run_phenotypic_mean',
# 	outputDir = '\\\\d\\dfs\\Groups\\biol\\sysbc\\claassen\\macnairw\\Outputs\\08_treeSNE\\Bodenmiller\\bodenmiller_test_run_phenotypic'
# 	)
# clusterParams = list(seed = 0, approach = 'kmeans', krange = 3:10, nstarts = 5000)
# # clusterParams = list(seed = 0, approach = 'kmeans', k = 15, nstarts = 5000)
# # clusterParams = list(seed = 0, approach = 'pam', krange = 1:40)
# treeSNE_clusterer(runParams, clusterParams)


runParams = list(
	saveStem = 'treeSNE_test_prostate_3_mean',
	outputDir = '/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/treeSNE_test_prostate_3'
	)
# clusterParams = list(seed = 0, approach = 'kmeans', krange = c(3:9), nstarts = 5000)
clusterParams = list(seed = 0, approach = 'kmeans', k = 8, nstarts = 5000)
# clusterParams = list(seed = 0, approach = 'pam', krange = 1:25)
treeSNE_clusterer(runParams, clusterParams)
