function [FlowSPD_output_filename] = for_4_Process_DataForSPADE_replace_GUI_treeSNE_randomised(treeSNE_parameters)

loaded_variables = load_variables(treeSNE_parameters);

loaded_variables = calc_clusters_and_MST(loaded_variables);

save_outputs(loaded_variables);

FlowSPD_output_filename = loaded_variables.result_filename;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load_variables:
function [loaded_variables] = load_variables(treeSNE_parameters)

	% load variables from previous function
	loaded_variables.data_filename = treeSNE_parameters.FlowSPADE_input_filename;
	load(loaded_variables.data_filename)
	
	% load data
	loaded_variables.data = data; 
	
	% check/load marker names
	if exist('marker_names') && length(marker_names)==size(data,1)
		loaded_variables.marker_names = marker_names;
	else
		error('marker_names does not exist, or its length does not match the data');
	end

	% check/load local density
	if exist('local_density')
		loaded_variables.local_density = local_density;
	else
		error('local_density does not exist');
	end

	% check/load used markers
	if exist('used_markers')
		loaded_variables.used_markers = used_markers;
	else
		error('used_markers does not exist');
	end

	% check/load extra markers
	if exist('extra_markers')
		loaded_variables.extra_markers = extra_markers;
	else
		error('extra_markers does not exist');
	end

	loaded_variables.target_num_sample_groups = treeSNE_parameters.num_target_clusters;
	loaded_variables.sample_group_assign = [];
	
	% set save_filename
    loaded_variables.result_filename = [treeSNE_parameters.runStem '_clusters.mat'];
	
	return


%% calc_clusters_and_MST:
function [loaded_variables] = calc_clusters_and_MST(loaded_variables)

	% grouping and find mst
	[mst_adj_dist, sample_group_assign, sample_modules] = ...
		group_flow_samples_to_tree_treeSNE_randomised(loaded_variables.data, loaded_variables.used_markers, loaded_variables.extra_markers, loaded_variables.target_num_sample_groups);
	
	% put outputs into variable
	loaded_variables.mst_adj_dist 			= mst_adj_dist;
	loaded_variables.sample_group_assign 	= sample_group_assign;
	loaded_variables.sample_modules 		= sample_modules;

	% remove unassigned groups
	loaded_variables.data(:,loaded_variables.sample_group_assign==0)				=[];
	loaded_variables.local_density(loaded_variables.sample_group_assign==0)			=[];
	loaded_variables.sample_group_assign(loaded_variables.sample_group_assign==0)	=[];

	return


%% save_outputs:
function [] = save_outputs(loaded_variables)
	% set up variables to be saved
	data 						= loaded_variables.data;
	marker_names 				= loaded_variables.marker_names;
	used_markers 				= loaded_variables.used_markers;
	extra_markers 				= loaded_variables.extra_markers;
	target_num_sample_groups 	= loaded_variables.target_num_sample_groups;
	sample_group_assign 		= loaded_variables.sample_group_assign;
	mst_adj_dist 				= loaded_variables.mst_adj_dist;
	data_filename 				= loaded_variables.data_filename;
	result_filename 			= loaded_variables.result_filename;
	local_density 				= loaded_variables.local_density;
	sample_modules 				= loaded_variables.sample_modules;

	% save
	save(result_filename, 'data', 'marker_names', 'used_markers', 'extra_markers', 'target_num_sample_groups', ...
		'sample_group_assign', 'mst_adj_dist',  'data_filename', 'result_filename', 'local_density', ...
		'sample_modules');

	return
