function [] = treeSNE_MST_section_randomised(FlowSPADE_output_filename, treeSNE_parameters, referenceCellIdx)
% treeSNE_tSNE_section.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDims 		= treeSNE_parameters.nDims;
saveStem 	= treeSNE_parameters.saveStem;
dataDir 	= treeSNE_parameters.dataDir;
outputDir 	= treeSNE_parameters.outputDir;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(outputDir);
[clusters sampleIdx local_density used_marker_names used_modules extra_marker_names extra_modules mst_adj_dist] = load_data(FlowSPADE_output_filename);

total_density = calc_total_density_by_sample(clusters, sampleIdx, local_density);

disp('Calculating MST distance...');
nn = size(mst_adj_dist,1);
mst_adj_dist = sparse(mst_adj_dist);
mst_dist = dijkstra(mst_adj_dist, (1:nn)');
% mst_dist = zeros(size(mst_adj_dist)); 			% this line is useful for code testing

% select reference clusters
referenceClusters = clusters(referenceCellIdx);

% find distances between cells under MST distances
referenceDists			= calc_reference_distances(referenceClusters, mst_dist);
referenceWeights		= calc_reference_weights(referenceClusters, total_density);
referenceUsedMarkers	= calc_reference_used_markers(referenceClusters, used_modules);
referenceExtraMarkers	= calc_reference_extra_markers(referenceClusters, extra_modules);

% save outputs
save_outputs(referenceWeights, referenceDists, referenceUsedMarkers, used_marker_names, referenceExtraMarkers, extra_marker_names, treeSNE_parameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load_data: function description
function [clusters sampleIdx local_density used_marker_names used_modules extra_marker_names extra_modules mst_adj_dist] = load_data(FlowSPADE_output_filename)
	load(FlowSPADE_output_filename);
	
	sampleIdx = data(end,:)';
	used_marker_names = marker_names(used_markers);
	extra_marker_names = marker_names(extra_markers);
	sample_modules = sample_modules';
	used_modules = sample_modules(:, used_markers);
	extra_modules = sample_modules(:, extra_markers);
	clusters = sample_group_assign';

	return


%% calc_total_density_by_sample: 
function [total_density] = calc_total_density_by_sample(clusters, sampleIdx, local_density)
	if numel(sampleIdx) ~= numel(clusters)
		error('Data and cluster vector lengths not equal.');
	end

	uniqueSamples = unique(sampleIdx); nSamples = numel(uniqueSamples);
	uniqueClusters = unique(clusters); nClusters = numel(uniqueClusters);

	total_density = zeros(nClusters, nSamples);
	% parfor ii = 1:nSamples
	for jj = 1:nClusters
		clusterIdx = (clusters == uniqueClusters(jj));
		for ii = 1:nSamples
			sample = (sampleIdx == uniqueSamples(ii));
			idx = and(sample,clusterIdx);
			total_density(jj, ii) = sum(local_density(idx));
		end % ii
	end % jj

	return


%% calc_reference_distances: calculates distances between all pairs of reference cells, using MST cluster distance as proxy
function [referenceDists] = calc_reference_distances(referenceClusters, mst_dist)
	referenceDists = mst_dist(referenceClusters, referenceClusters);

%% calc_reference_distances: calculates distances between all pairs of reference cells, using MST cluster distance as proxy
function [referenceWeights] = calc_reference_weights(referenceClusters, total_density);
	referenceWeights = total_density(referenceClusters,:);

%% calc_reference_distances: calculates distances between all pairs of reference cells, using MST cluster distance as proxy
function [referenceUsedMarkers] = calc_reference_used_markers(referenceClusters, used_modules);
	referenceUsedMarkers = used_modules(referenceClusters,:);

%% calc_reference_distances: calculates distances between all pairs of reference cells, using MST cluster distance as proxy
function [referenceExtraMarkers] = calc_reference_extra_markers(referenceClusters, extra_modules);
	referenceExtraMarkers = extra_modules(referenceClusters,:);


%% save_outputs: 
function [] = save_outputs(referenceWeights, referenceDists, referenceUsedMarkers, used_marker_names, referenceExtraMarkers, extra_marker_names, treeSNE_parameters)
	savePath = fullfile(treeSNE_parameters.outputDir, [treeSNE_parameters.runStem '_outputs.mat']);
	save(savePath, 'referenceWeights', 'referenceDists', 'referenceUsedMarkers', 'referenceExtraMarkers');

% %% save_txt_file: 
% function [] = save_txt_file(saveFilename, header, saveData)
	
% 	if isempty(header)

% 		disp(['Saving file ' saveFilename]);

% 		dataSpec = ['%4.4f' repmat('\t%4.4f', 1, size(saveData,2) -1) '\n'];;
% 		fid = fopen(saveFilename, 'w');
% 		for ii = 1:size(saveData,1)
% 			fprintf(fid, dataSpec, saveData(ii,:));
% 		end
% 		fclose(fid);

% 	else

% 		if size(header,2) ~= size(saveData,2)
% 			error(['Problem saving ' saveFilename ': header and data are not compatible lengths.']);
% 		else
% 			disp(['Saving file ' saveFilename]);
% 		end

% 		headerSpec = ['%s' repmat('\t%s', 1, size(header,2) -1) '\n'];
% 		dataSpec = ['%4.4f' repmat('\t%4.4f', 1, size(saveData,2) -1) '\n'];;
% 		fid = fopen(saveFilename, 'w');
% 		fprintf(fid, headerSpec, header{:});
% 		for ii = 1:size(saveData,1)
% 			fprintf(fid, dataSpec, saveData(ii,:));
% 		end
% 		fclose(fid);

% 	end % if


% %% save_outputs: 
% function [] = save_outputs(referenceWeights, referenceDists, referenceUsedMarkers, used_marker_names, referenceExtraMarkers, extra_marker_names, treeSNE_parameters)
% % function [] = save_outputs(referenceMarkers, marker_names, treeSNE_parameters)

% 	cd(treeSNE_parameters.outputDir);
% 	runStem = treeSNE_parameters.runStem;

% 	% set up density save
% 	saveFilename = [runStem '_weights' '.txt'];
% 	header = treeSNE_parameters.file_annot';
% 	saveData = referenceWeights;
% 	% save density
% 	save_txt_file(saveFilename, header, saveData);

% 	% set up distance save
% 	saveFilename = [runStem '_MST_dists' '.txt'];
% 	header = [];
% 	saveData = referenceDists;
% 	% save distance
% 	save_txt_file(saveFilename, [], saveData);

% 	% set up module save
% 	saveFilename = [runStem '_used_markers' '.txt'];
% 	header = used_marker_names';
% 	saveData = referenceUsedMarkers;
% 	% save used modules
% 	save_txt_file(saveFilename, header, saveData);

% 	% set up module save
% 	saveFilename = [runStem '_extra_markers' '.txt'];
% 	header = extra_marker_names';
% 	saveData = referenceExtraMarkers;
% 	% save extra modules
% 	save_txt_file(saveFilename, header, saveData);

