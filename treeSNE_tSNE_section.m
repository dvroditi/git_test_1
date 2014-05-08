function [] = treeSNE_tSNE_section(FlowSPADE_output_filename, treeSNE_parameters)
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
nn=size(mst_adj_dist,1);
mst_dist = dijkstra(mst_adj_dist, (1:nn)');
% mst_dist = zeros(size(mst_adj_dist)); 			% this line is useful for code testing
disp('Calculating tSNE...');
tsne_outputs = tsne_d(mst_dist, [], 2);
% tsne_outputs = zeros(size(mst_dist,1), nDims);	% this line is useful for code testing

save_outputs(total_density, tsne_outputs, used_modules, used_marker_names, extra_modules, extra_marker_names, treeSNE_parameters);
disp(['treeSNE completed for ' saveStem]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%q

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

	return


%% save_outputs: 
function [] = save_outputs(total_density, tsne_outputs, used_modules, used_marker_names, extra_modules, extra_marker_names, treeSNE_parameters)

	saveStem = treeSNE_parameters.saveStem;

	% set up density save
	saveFilename = [saveStem '_weights' '.txt'];
	header = treeSNE_parameters.file_annot';
	saveData = total_density;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up tsne save
	saveFilename = [saveStem '_tsne' '.txt'];
	header = concatenate_cell_strings('Y', (1:treeSNE_parameters.nDims)');
	saveData = tsne_outputs;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up used module save
	saveFilename = [saveStem '_used_markers' '.txt'];
	header = used_marker_names';
	saveData = used_modules;
	% save density
	save_txt_file(saveFilename, header, saveData);

	% set up used module save
	saveFilename = [saveStem '_extra_markers' '.txt'];
	header = extra_marker_names';
	saveData = extra_modules;
	% save density
	save_txt_file(saveFilename, header, saveData);

	return

