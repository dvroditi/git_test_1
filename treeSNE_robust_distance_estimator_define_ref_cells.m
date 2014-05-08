
% Reference cell index needed; other values obsolete 
% Previously information on position of ref cell in marker space was purely based on ref cells;
% Now ref cell position in marker space is represented by mean marker value across cluster centres, for the clusters containing the ref cell (to be used in classifier)

% Saves information in marker space about reference cells 
function [] = treeSNE_robust_distance_estimator_define_ref_cells(treeSNE_parameters)

	% set seeds
	stream0 = RandStream('mt19937ar', 'Seed', treeSNE_parameters.refSeed);
	RandStream.setGlobalStream(stream0);

	% set up variables
	outputDir = treeSNE_parameters.outputDir;
	cd(outputDir);
	
	% load mat file
	FlowSPD_input_filename = 'DataForSPADE.mat';
	load(FlowSPD_input_filename);

	% take sample
	referenceCells = sort(randsample(size(data,2), treeSNE_parameters.num_sample_cells));
	referenceCellIdx = ismember(1:size(data,2), referenceCells);

	% need: used_markers (as marker_names); cell values
	used_marker_names = marker_names(used_markers);
	used_modules = data(used_markers, referenceCellIdx)';
	
	extra_marker_names = marker_names(extra_markers);
	extra_modules = data(extra_markers, referenceCellIdx)';
	
	% save
	referenceFile = [treeSNE_parameters.saveStem '_reference_cells.mat'];
	save(referenceFile, 'referenceCellIdx', 'used_marker_names', 'used_modules', 'extra_marker_names', 'extra_modules');
