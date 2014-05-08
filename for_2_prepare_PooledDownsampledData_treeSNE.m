function [] = for_2_prepare_PooledDownsampledData_treeSNE(treeSNE_parameters)

% set up input variables
[filenames target_prctile exclude_prctile PooledDownsampledDataFilename] = define_parameters(treeSNE_parameters);

% set up output variables
all_data=[];
file_indices=[];
all_local_density=[];
used_marker_names = [];
extra_marker_names = [];
all_marker_names = [];
RefDataSize = [];

% loop through files
for ii=1:length(filenames)

	% load file
	say_loading_file(ii, filenames);
	load([filenames{ii}(1:end-4),'.mat'])

	% update used_marker_names
	[RefDataSize, all_marker_names, used_marker_names, extra_marker_names] = update_marker_names(ii, RefDataSize, data, all_marker_names, marker_names, used_marker_names, used_markers, extra_marker_names, extra_markers);
	
	% exclude low density points
	[data local_density] = exclude_low_density(data, local_density, exclude_prctile);

	% calculate target density
	target_density = calc_target_density(target_prctile, local_density);

	% do downsampling and renormalizing local density
	[data local_density] = downsample_file(target_density, local_density, data, RefDataSize);

	% store values, mess around with marker names
	[all_data all_marker_names all_local_density file_indices] = add_new_values_to_storage_variables(data, marker_names, local_density, all_data, all_marker_names, all_local_density, file_indices, ii);

end

% tidy up
[data local_density marker_names used_markers extra_markers] = tidy_up(all_data, file_indices, all_marker_names, all_local_density, marker_names, used_markers, used_marker_names, extra_markers, extra_marker_names);

% save
say_saving_file(data, filenames);
save(PooledDownsampledDataFilename, 'data', 'local_density', 'marker_names', 'used_markers', 'extra_markers');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define_parameters:
function [filenames target_prctile exclude_prctile PooledDownsampledDataFilename] = define_parameters(treeSNE_parameters)
	
	filenames = treeSNE_parameters.filenames;
	target_prctile = treeSNE_parameters.target_prctile;
	exclude_prctile = treeSNE_parameters.exclude_prctile;
	PooledDownsampledDataFilename = treeSNE_parameters.PooledDownsampledDataFilename;
	
	return

%% say_loading_file:
function [] = say_loading_file(ii, filenames)
	display(['downsampling and pooling fcs file: ',num2str(ii),'/',num2str(length(filenames))]);
	display(filenames{ii});
	
	return


%% update_marker_names: 
function [RefDataSize, all_marker_names, used_marker_names, extra_marker_names] = update_marker_names(ii, RefDataSize, data, all_marker_names, marker_names, used_marker_names, used_markers, extra_marker_names, extra_markers)

	if ii==1
		RefDataSize = size(data,2);
		all_marker_names = marker_names;
		used_marker_names = marker_names(used_markers);
		extra_marker_names = marker_names(extra_markers);
	else
		RefDataSize = RefDataSize;
		all_marker_names = all_marker_names;
		used_marker_names = unique([used_marker_names; marker_names(used_markers)]);
		extra_marker_names = unique([extra_marker_names; marker_names(extra_markers)]);
	end

	return


%% exclude_low_density:
function [data local_density] = exclude_low_density(data, local_density, exclude_prctile)

	% SPADE does this with <=, but I found this gave odd behaviour (e.g. entirely excluding a file)
	data(:,local_density<prctile(local_density,exclude_prctile))=[];
	local_density(local_density<prctile(local_density,exclude_prctile))=[];

	return


%% calc_target_density:
function [target_density] = calc_target_density(target_prctile, local_density)
	
	if target_prctile<100
		target_density = prctile(local_density,target_prctile);
	elseif target_prctile>=100  % then this variable contains the number of desired cells we want after downsampling
		num_desired_cells = target_prctile;
		target_density = calc_target_density_for_desired_num_of_cells(local_density, num_desired_cells);
	end

	return


%% calc_target_density_for_desired_num_of_cells:
function [target_density] = calc_target_density_for_desired_num_of_cells(local_density, desired_num)
	% keep_prob = x./local_density
	% need to find the value of "x", such that if we downsample according to
	% "keep_prob", we end up with about "desired_num" cells
	% therefore, need to solve the following
	%      sum(min(x/local_density(i),1)) = desired_num
	% which is equivalent to
	%      x = (desired_num-i) / sum(1/local_density(i+1:end)) && local_density(i)<=x<=local_density(i+1) 

	if desired_num>=length(local_density)
		target_density = max(local_density)+1;
		return
	end
	ld = [sort(local_density,'ascend')];
	if desired_num/sum(1./local_density) <= ld(1)
		target_density = desired_num/sum(1./local_density);
		return
	end
	for i=1:length(ld)-1
		x = (desired_num-i) / sum(1./ld(i+1:end));
		if ld(i)<=x && x<=ld(i+1) 
			break;
		end
	end
	target_density = x;
	return


%% downsample_file: 
function [data local_density] = downsample_file(target_density, local_density, data, RefDataSize)

	% decide which cells have to go
	keep_prob = min(1,(target_density./local_density));
	is_keep = rand(1,length(local_density))<keep_prob;  
	% remove cells with dodgy values
	is_keep(find(sum(isnan(data))~=0))=0;
	display([num2str(sum(is_keep)),' cells kept in this fcs file'])
	% exclude other cells
	data = data(:,is_keep);
	% renormalize local_density
	local_density = local_density(is_keep)/length(is_keep)*RefDataSize;

	return

%% add_new_values_to_storage_variables:
function [all_data all_marker_names all_local_density file_indices] = add_new_values_to_storage_variables(data, marker_names, local_density, all_data, all_marker_names, all_local_density, file_indices, ii)
	if isequal(marker_names,all_marker_names)
		all_data = [all_data,data];
	else
		new_marker_names = setdiff(marker_names,all_marker_names);
		all_marker_names = [all_marker_names;new_marker_names];
		all_data = [all_data;repmat(NaN,length(new_marker_names),size(all_data,2))];
		data_tmp = zeros(size(all_data,1),size(data,2))+NaN;
		[C,IA,IB] = intersect(marker_names,all_marker_names);
		data_tmp(IB,:) = data(IA,:);
		all_data = [all_data, data_tmp];
	end

	all_local_density = [all_local_density,local_density];
	file_indices = [file_indices,repmat(ii,1,size(data,2))];

	return

%% tidy_up:
function [data local_density marker_names used_markers extra_markers] = tidy_up(all_data, file_indices, all_marker_names, all_local_density, marker_names, used_markers, used_marker_names, extra_markers, extra_marker_names)
	% add file indices to data matrix
	data = [all_data;file_indices];
	all_marker_names{end+1} = 'FileInd';

	marker_names = all_marker_names;
	local_density = all_local_density;

	[C,used_markers,IB] = intersect(marker_names, used_marker_names); 
	used_markers = sort(used_markers);

	[C,extra_markers,IB] = intersect(marker_names, extra_marker_names); 
	extra_markers = sort(extra_markers);

	return


%% say_saving_file:
function [] =say_saving_file(data, filenames);
	display(' ')
	display(['PooledDownsampledData has ', num2str(size(data,2)), ' cells from ', num2str(length(filenames)), ' files'])
	display(' ')

	return

