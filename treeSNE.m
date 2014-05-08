function [] = treeSNE(treeSNE_parameters)

	addpath('/Users/delaura/Documents/MATLAB/Code/treeSNE');
	addpath('/Users/delaura/Documents/DATA/PC_CyTOF/140212_AJa57/experiment_958_files');
		

	% check parameters are ok
	check_treeSNE_parameters(treeSNE_parameters);

	% SPADE parameters
	treeSNE_parameters.PooledDownsampledDataFilename	= 'PooledData.mat';
	treeSNE_parameters.FlowSPADE_input_filename			= 'DataForSPADE.mat';
	outputDir = treeSNE_parameters.outputDir;

	% get filenames we will use
	cd(treeSNE_parameters.dataDir);
	filenames = treeSNE_parameters.filenames;

	% for each fcs file, compute the local density of each cell
	% save the data and the local density in an intermediate .mat file
	% one intermediate .mat file is created for each fcs file
	% (20 min per file if use parallel computing toolbox; if not, ~60 mins)
	if treeSNE_parameters.pool_flag
		if ~matlabpool('size')
			matlabpool;
		end
		parfor i=1:length(filenames)
			filename = filenames{i};
			if exist(fullfile(outputDir, [filenames{i}(1:end-4),'.mat']))==2
				continue;
			end
			tic 
			parfor_1_read_individual_fcs_get_downsampling_info_treeSNE(filename, treeSNE_parameters); 
			toc
		end % for
	else
		for i=1:length(filenames)
			filename = filenames{i};
			if exist(fullfile(outputDir, [filenames{i}(1:end-4),'.mat']))==2
				continue;
			end
			tic 
			for_1_read_individual_fcs_get_downsampling_info_treeSNE(filename, treeSNE_parameters);
			toc
		end % for
	end % if

	cd(outputDir);
	% This script actually does the downsampling. Each fcs file is downsampled separately, 
	% all the downsampled cells are pooled into a single file, generating "PooledData.mat"
	% although named as parfor_..., this function does not need parallel computing toolbox
	for_2_prepare_PooledDownsampledData_treeSNE(treeSNE_parameters);

	% prepare input file for the clustering step
	% if the total number of pooled cells is too large, we uniformly downsample the data 
	% further, because the current clustering algorithm cannot handle too many cells without 
	% memory problem or running forever. We limit the total number to be 50000.
	for_3_Prepare_DataForSPADE_treeSNE(treeSNE_parameters);


	% clustering and MST construction
	% this step generate a .mat file that has a very long name
	[FlowSPADE_output_filename] = for_4_Process_DataForSPADE_replace_GUI_treeSNE(treeSNE_parameters);

	% do mst-based tsne
	treeSNE_tSNE_section(FlowSPADE_output_filename, treeSNE_parameters);
