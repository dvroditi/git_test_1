% dijkstra is a slow step that can be sped up using C - called within treeSNE_MST_section_randomised
% Calculates random clustering of data and MST between clusters and records for each ref cell: 
% - location, size and composition of the cluster it was put in
% - distances to other ref cells, by going along MST
function [] = treeSNE_robust_distance_estimator(treeSNE_parameters, seed) 

	% set seed
	stream0 = RandStream('mt19937ar','Seed',seed);
	RandStream.setGlobalStream(stream0);

	% set up variables
	treeSNE_parameters.runStem = [treeSNE_parameters.saveStem '_cluster_run_' sprintf('%04d', seed)];
	cd(treeSNE_parameters.outputDir);
	
	% do randomised clustering
	% FlowSPADE_output_filename = [treeSNE_parameters.runStem '_clusters.mat'];
	[FlowSPADE_output_filename] = for_4_Process_DataForSPADE_replace_GUI_treeSNE_randomised(treeSNE_parameters);
	
	% load reference cell index
	referenceFile = [treeSNE_parameters.saveStem '_reference_cells.mat'];
	load(referenceFile);
	
	% do mst
	treeSNE_MST_section_randomised(FlowSPADE_output_filename, treeSNE_parameters, referenceCellIdx);

