% mst_dist = dijkstra(mean_distance, []); Slow step that can be sped up using C
% Averaging of the clusters where ref cells fell into - includes: location, size, composition and intercluster (tree-induced) distances
% Runs T-SNE
function [] = treeSNE_robust_distance_aggregator(treeSNE_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% set up variables
	nDims = treeSNE_parameters.nDims;
	saveStem = treeSNE_parameters.saveStem;
	outputDir = treeSNE_parameters.outputDir;
	num_sample_cells = treeSNE_parameters.num_sample_cells;
	nRuns = treeSNE_parameters.nRuns;

	nSamples = numel(treeSNE_parameters.file_annot);
	nUsedMarkers = numel(treeSNE_parameters.used_markers);
	nExtraMarkers = numel(treeSNE_parameters.extra_markers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADMIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	cd(outputDir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% how to store? reshape into columns
	all_distances = zeros(num_sample_cells.^2, nRuns);
	all_weights = zeros(num_sample_cells*nSamples, nRuns);
	all_used_locations = zeros(num_sample_cells*nUsedMarkers, nRuns);
	all_extra_locations = zeros(num_sample_cells*nExtraMarkers, nRuns);

	% open all files
	disp('Loading all runs to calculate mean values...');
	% idx = setdiff(1:nRuns, [74, 83, 85, 87]);
	idx = 1:nRuns;
	for ii = idx
		
		runStem = [saveStem '_cluster_run_' sprintf('%04d', ii)];
		
		[run_distances run_weights run_used_locations run_extra_locations] = get_data(runStem, treeSNE_parameters);

		all_distances(:, ii) 		= reshaper(run_distances);
		all_weights(:, ii) 			= reshaper(run_weights);
		all_used_locations(:, ii) 	= reshaper(run_used_locations);
		all_extra_locations(:, ii) 	= reshaper(run_extra_locations);

	end

	% calculate and plot stds
	plot_stds(all_distances, all_weights);
	
	% calculate means
	mean_distance 		= reshape(mean(all_distances, 2), num_sample_cells, num_sample_cells);
	mean_weight 			= reshape(mean(all_weights, 2), num_sample_cells, nSamples);
	mean_used_location	= reshape(mean(all_used_locations, 2), num_sample_cells, nUsedMarkers);
	mean_extra_location	= reshape(mean(all_extra_locations, 2), num_sample_cells, nExtraMarkers);

	% do mst, distance etc
	disp('Calculating mst distance for mean distances...');
    
    nn=size(mst_adj_dist,1);
    
% 	mst_dist = dijkstra(mean_distance, []); % makes it a full matrix of distances between all cells
	mst_dist = dijkstra(mean_distance,(1:nn)');

    disp('Calculating tSNE for mean distances...');
	tsne_outputs = tsne_d(mst_dist, [], 2);
	% load other inputs
	load([saveStem '_reference_cells.mat']);

	save_outputs(mean_weight, tsne_outputs, mean_used_location, used_marker_names, mean_extra_location, extra_marker_names, treeSNE_parameters);
	disp(['Robustly estimated treeSNE completed for ' saveStem]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get_data: 
function [output] = reshaper(inputthing)
	output = reshape(inputthing, size(inputthing,1)*size(inputthing,2), 1);


%% get_data: 
function [run_distances run_weights run_used_locations run_extra_locations] = get_data(runStem, treeSNE_parameters)
	
	filepath = fullfile(treeSNE_parameters.outputDir, [runStem '_outputs.mat']);
	load(filepath);

	run_distances = referenceDists;
	run_weights = referenceWeights;
	run_used_locations = referenceUsedMarkers;
	run_extra_locations = referenceExtraMarkers;

	return

%% plot_stds: plots histograms of std
function [] = plot_stds(all_distances, all_weights)
	distanceStd = std(all_distances, 0, 2);
	subplot(1,2,1), hist(distanceStd);

	weightStd = std(all_weights, 0, 2);
	subplot(1,2,2), hist(weightStd);

%% save_txt_file: 
function [] = save_txt_file(saveFilename, header, saveData)
	if size(header,2) ~= size(saveData,2)
		error(['Problem saving ' saveFilename ': header and data are not compatible lengths.']);
	else
		disp(['Saving file ' saveFilename]);
	end

	headerSpec = ['%s' repmat('\t%s', 1, size(header,2) -1) '\n'];
	dataSpec = ['%4.4f' repmat('\t%4.4f', 1, size(saveData,2) -1) '\n'];;
	fid = fopen(saveFilename, 'w');
	fprintf(fid, headerSpec, header{:});
	for ii = 1:size(saveData,1)
		fprintf(fid, dataSpec, saveData(ii,:));
	end
	fclose(fid);

%% save_outputs: 
function [] = save_outputs(mean_weight, tsne_outputs, mean_used_location, used_marker_names, mean_extra_location, extra_marker_names, treeSNE_parameters)
% function [] = save_outputs(meanLocation, marker_names, treeSNE_parameters)

	cd(treeSNE_parameters.outputDir);
	saveStem = treeSNE_parameters.saveStem;

	% set up density save
	saveFilename = [saveStem '_mean_weights' '.txt'];
	header = treeSNE_parameters.file_annot';
	saveData = mean_weight;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up tsne save
	saveFilename = [saveStem '_mean_tsne' '.txt'];
	header = concatenate_cell_strings('Y', (1:treeSNE_parameters.nDims)');
	saveData = tsne_outputs;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up module save
	saveFilename = [saveStem '_mean_used_markers' '.txt'];
	header = used_marker_names';
	saveData = mean_used_location;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up module save
	saveFilename = [saveStem '_mean_extra_markers' '.txt'];
	header = extra_marker_names';
	saveData = mean_extra_location;
	% save density
	save_txt_file(saveFilename, header, saveData);
