function for_1_read_individual_fcs_get_downsampling_info_treeSNE(filename, treeSNE_parameters)

	[data fcshdr marker_names used_markers extra_markers] = import_data(filename, treeSNE_parameters);

	[local_density kernel_width] = estimate_local_density(data(used_markers,:), treeSNE_parameters.kernel_width_para);

	filePath = fullfile(treeSNE_parameters.outputDir, [filename(1:end-4),'.mat']);
	save(filePath, 'data', 'marker_names', 'local_density', 'used_markers', 'extra_markers', 'kernel_width');

	return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% import_data: 
function [data fcshdr marker_names used_markers extra_markers] = import_data(filename, treeSNE_parameters)

	used_markers = treeSNE_parameters.used_markers;
	extra_markers = treeSNE_parameters.extra_markers;

	[fcsdat, fcshdr, fcsdatscaled] = fca_readfcs_treeSNE(filename, treeSNE_parameters.fcs_cell_limit);
	data = flow_arcsinh(fcsdat', treeSNE_parameters.arcsinh_cofactor);
	data = data(:,1:min(end,500000));

	display(['Read and downsample fcs file: '])
	display(['   ',filename,]);
	display(['   ',num2str(size(data,2)), ' cells']);
	% get marker names from file
	for i=1:length(fcshdr.par), 
		marker_names{i,1} = fcshdr.par(i).name2;  
		if isequal(unique(fcshdr.par(i).name2),' ') || isempty(fcshdr.par(i).name2)
			marker_names{i,1} = fcshdr.par(i).name;
		end
	end
	
	% normalize
	if exist('is_normalize') && is_normalize==1
		if ~exist('normalize_weight_factor')
			data = flow_normalize(data);
		else
			data = flow_normalize(data, normalize_weight_factor);
		end
	end
	
	% find index of used markers
	if iscell(used_markers)
		[C,IA,IB] = intersect(used_markers,marker_names);
		if length(used_markers)~=length(C)
			error('Some used markers do not exist in at least one FCS file');
		end
		used_markers = sort(IB);
	end

	% find index of extra markers
	if iscell(extra_markers)
		[C,IA,IB] = intersect(extra_markers,marker_names);
		if length(extra_markers)~=length(C)
			error('Some used markers do not exist in at least one FCS file');
		end
		extra_markers = sort(IB);
	end

	return


%% flow_arcsinh: calculates arcsinh values of data
function [data] = flow_arcsinh(data, cofactor)
	% cofactor can be one signal number, for all channels
	%          can be a vector, one element for one channel
	% if cofactor == 0, then it means linear, no transform. 

	if length(cofactor)~=1 && length(cofactor)~=size(data,1)
		disp('wrong input of cofactors')
		disp('data not transformed')
		return
	end

	if length(cofactor)==1
		if cofactor==0
			return
		else
			data = log( data(:,:)/cofactor + sqrt((data(:,:)/cofactor).^2+1) );
			return
		end
	end

	if length(cofactor) == size(data,1)
		for i=1:size(data,1)
			if cofactor(i)==0
				continue;
			else
				data(i,:) = log( data(i,:)/cofactor(i) + sqrt((data(i,:)/cofactor(i)).^2+1) );
			end        
		end
	end

	return


%% flow_normalize: normalizes data so that 0 = 5th percentile, 1 = 95th percentile
function data = flow_normalize(data, normalize_weight_factor)

for i=1:size(data,1)
    [ll] = prctile(data(i,:),5);
    [uu] = prctile(data(i,:),95);
    data(i,:) = (data(i,:)-ll)./(uu-ll)*5;
end

if exist('normalize_weight_factor')
    for i=1:size(data,1)
        data(i,:) = data(i,:) * normalize_weight_factor(i);
    end
end
return


%% estimate_local_density: downsamples this file
function [local_density kernel_width] = estimate_local_density(data, kernel_width_para)

	% if parameters not already defined, use defaults
	if ~exist('kernel_width_para'), kernel_width_para = 10; end
	% if ~exist('target_prctile'),    target_prctile = 5;     end

	% estimate median distance to nearest neighbour
	med_min_dist = calc_med_min_dist(data);

	% define kernel width with respect to minimum distance
	kernel_width = kernel_width_para * med_min_dist;
	fprintf('   For this %d channel data, KERNEL_WIDTH is %3.3f\n', size(data,1),kernel_width);

	[local_density] = calc_local_density(data, kernel_width, med_min_dist);

	% this section calculates keeping probabilities, but they not used in SPADE
	% target_density = prctile(local_density,target_prctile);
	% keep_prob = (target_density./local_density);
	% fprintf('\n   Down-sample cells ... \n');
	% is_keep = rand(1,size(data,2))<keep_prob;

	return


%% calc_med_min_dist: calculates distance to nearest point for sample of 2000 points, returns median
function [med_min_dist] = calc_med_min_dist(data)

	fprintf('   finding empirical dist of the min distance between cells ... \n');
	min_dist = zeros(1,min([size(data,2),2000,floor(2500000000/size(data,2))]));  % (1) number of cells if data matrix is small, (2) block size = 2000 if data large, (3) if number of cells is extremely large, we need to limit block size, otherwise line 38 will be out of memory

	ind = sort(randsample(1:size(data,2),length(min_dist)));
	data_tmp = data(:,ind);
	all_dist = zeros(length(ind),size(data,2));
	for i=1:size(data,2)
		all_dist(:,i) = sum(abs(repmat(data(:,i),1,size(data_tmp,2)) - data_tmp),1)';
	end
	all_dist(all_dist==0)=max(max(all_dist));

	min_dist = min(all_dist');

	med_min_dist = median(min_dist);

	return

%% calc_local_density: calculates local density of data
function [local_density] = calc_local_density(data, kernel_width, med_min_dist)

	local_density = zeros(1,size(data,2));
	local_density_tmp = zeros(1,size(data,2));
	fprintf('   finding local density for each cell ... \n   %10d',0);
	count = 1; 
	while sum(local_density==0)~=0
		ind = find(local_density==0); ind = ind(1:min(1000,end));
		data_tmp = data(:,ind);
		local_density_tmp = local_density(ind);
		all_dist = zeros(length(ind),size(data,2));
		for i=1:size(data,2)
			all_dist(:,i) = sum(abs(repmat(data(:,i),1,size(data_tmp,2)) - data_tmp),1)';
		end
		for i=1:size(data_tmp,2)
			% local_density_tmp(i) = sum(kernel_width>=all_dist(i,:)); % sum(max(kernel_width-all_dist(i,:),0)/kernel_width); % 
			% local_density(all_dist(i,:) < 1.5*med_min_dist) = local_density_tmp(i);
			local_density_tmp(i) = sum(max(kernel_width-all_dist(i,:),0)/kernel_width);
			local_density(all_dist(i,:) < med_min_dist) = local_density_tmp(i);
		end
		fprintf('\b\b\b\b\b\b\b\b\b\b%10d',sum(local_density~=0));
	end

	return