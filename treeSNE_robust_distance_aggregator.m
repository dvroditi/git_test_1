function [] = treeSNE_robust_distance_aggregator(treeSNE_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    addpath('/Users/delaura/Documents/MATLAB/Code/treeSNE');
% 	addpath('/Users/delaura/Documents/DATA/PC_CyTOF/140212_AJa57/experiment_958_files');

    addpath('/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/prostate_run_split_large');	% location of data generated after first treeSNE_runs run
	
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
	% idx = [5 6 11 12 17 18 23 24 29 30 35 36 41 42 46 47 51 52 55 56 59 60 63 64];
	% idx = setdiff(1:nRuns, [5 6 11 12 17 18 23 24 29 30 34 35 36 41 42 45 46 47 51 52 55 56 59 60 63 64]);
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
	plot_stds(all_distances, all_weights, treeSNE_parameters);
	
	% calculate means
	mean_distance 		= reshape(mean(all_distances, 2), num_sample_cells, num_sample_cells);
	mean_weight 		= reshape(mean(all_weights, 2), num_sample_cells, nSamples);
	mean_used_location	= reshape(mean(all_used_locations, 2), num_sample_cells, nUsedMarkers);
	mean_extra_location	= reshape(mean(all_extra_locations, 2), num_sample_cells, nExtraMarkers);

	% do mst, distance etc
	% disp('Calculating mst distance for mean distances...');
	% mst_dist = dijkstra(mean_distance, []);
	% disp('Calculating tSNE for mean distances...');
	% tsne_outputs = tsne_d(mst_dist, [], 2);
	disp('Calculating tSNE for mean distances...');
	tsne_outputs = tsne_d(mean_distance, [], 2);
	% load other inputs
	load([saveStem '_reference_cells.mat']);

	save_outputs(mean_weight, tsne_outputs, mean_used_location, used_marker_names, mean_extra_location, extra_marker_names, treeSNE_parameters);
	disp(['Robustly estimated treeSNE completed for ' saveStem]);

	if treeSNE_parameters.delete_intermediate_files_when_finished
		delete_intermediate_files(treeSNE_parameters);
	end


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
function [] = plot_stds(all_distances, all_weights, treeSNE_parameters)
	distanceStd = std(all_distances, 0, 2);
	subplot(1,2,1), hist(distanceStd);

	weightStd = std(all_weights, 0, 2);
	subplot(1,2,2), hist(weightStd);

	saveas(gcf, [treeSNE_parameters.saveStem '_std_devs_of_runs.fig'], 'fig');
	

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

%% delete_intermediate_files: 
function [] = delete_intermediate_files(treeSNE_parameters)
	% set up variables
	n_runs = treeSNE_parameters.nRuns;
	file_list = {};

	% generate clusters filenames
	for ii = 1:n_runs
		filename = {[treeSNE_parameters.saveStem '_cluster_run_' sprintf('%04d', ii) '_clusters.mat']};
		file_list = [file_list; filename];
	end

	% generate outputs filenames
	for ii = 1:n_runs
		filename = {[treeSNE_parameters.saveStem '_cluster_run_' sprintf('%04d', ii) '_outputs.mat']};
		file_list = [file_list; filename];
	end

	% delete all filenames
	for filename = file_list
		delete(filename{:});
	end