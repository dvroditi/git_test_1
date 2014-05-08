function [mst_adj_dist, idx_samples, sample_modules] = group_flow_samples_to_tree_treeSNE_randomised(data, used_markers, extra_markers, target_num_sample_groups)

% randomise order of data
randOrder = randperm(size(data,2));
data = data(:, randOrder);


used_data = data(used_markers, :);

idx_samples = 1:size(data,2);
counter = 1;
while max(idx_samples)>=1.5*target_num_sample_groups
	fold = max(2,round(max(idx_samples)/5000)); 
	fprintf('merging %d groups into about %d groups, round %d\n', max(idx_samples), ceil(max(idx_samples)/fold), counter);
	tmp_ind = find(idx_samples~=0);
	if isempty(find(idx_samples==0))
		idx_samples = merge_sample_groups(used_data, fold, idx_samples);
	else
		idx_samples(tmp_ind) = merge_sample_groups(used_data(:,tmp_ind), fold,idx_samples(tmp_ind));
	end

% This section is in SPADE, and makes things faster by only doing 5 rounds then throwing remaining singletons away.
% This is unhelpful for making the outputs from multiple runs the same size, so I removed it.

%    if counter==5 % after 5 rounds of iterations, if a node still contains only one cell, this node is thrown away
%        for i=1:max(idx_samples)
%            if sum(idx_samples==i)==1
%                idx_samples(idx_samples==i)=0;
%            end
%        end
%        idx_samples = standardize_idx(idx_samples);
%        display(['throw away ',num2str(sum(idx_samples==0)),' singletons'])
%   end

	counter = counter + 1;
end 

% sample_modules = get_module_mean(data',idx_samples)';
sample_modules = calc_module_median(data, idx_samples);
used_sample_modules = sample_modules(used_markers, :);
[mst_adj mst_adj_dist] = mst(used_sample_modules');

% put data back into original order
idx_samples(randOrder) = idx_samples;

return



function idx_new = merge_sample_groups(data, fold, idx)
	% the samples grouped into give N groups by idx, 
	% reduce the number of groups by a factor of fold (min of fold is 2)
	% look at the groups one by one, 
	% for each group, use single linkage to find the closet other group
	% merge them together, and this other group is out of the game

	isactive_module = ones(1,max(idx));
	module_size = zeros(1,max(idx)); for i=1:length(module_size), module_size(i) = sum(idx==i); end
	isactive_sample = ones(1,size(data,2));
	fprintf('merging in progress ... groups left: %7d',max(idx));
	iter = 1; total_num_modules = max(idx);
	while sum(isactive_module)>1
		tmp = module_size; tmp(module_size==0)=Inf;  tmp(isactive_module==0) = Inf;
		candidate_module_ind = find(tmp==min(tmp(isactive_module==1)) & isactive_module==1);
		if length(candidate_module_ind)==1
			module_ind = candidate_module_ind;
		else
			module_ind = randsample(candidate_module_ind,1);
		end

		sample_in_module = find(idx==module_ind);
		module_center = median(data(:,sample_in_module),2);
		dist = sum(abs(repmat(module_center,1,size(data,2)) - data),1);
		dist(idx==module_ind) = Inf;    % dist to everyone in my own group
		[Y,I] = sort(dist,'ascend'); 
		module_to_be_deleted = idx(I(1:fold-1));
		first_inactive_module_ind = find(isactive_module(module_to_be_deleted)==0,1);
		if ~isempty(first_inactive_module_ind) && first_inactive_module_ind==1, module_to_be_deleted=[]; end
		if ~isempty(first_inactive_module_ind) && first_inactive_module_ind~=1, module_to_be_deleted=module_to_be_deleted(1:first_inactive_module_ind-1); end
		if isempty(module_to_be_deleted)
			isactive_module(module_ind)=0;
			fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d ',total_num_modules,sum(isactive_module==1));drawnow;
			continue;
		end
		isactive_module(module_to_be_deleted)=0;
		isactive_module(module_ind)=0;
		for i=1:length(module_to_be_deleted)
			idx(idx==module_to_be_deleted(i)) = module_ind;  % merge
		end
		module_size(module_to_be_deleted)=0; module_size(module_ind) = sum(idx==module_ind);
		isactive_sample(idx==module_ind)=0;
		total_num_modules = total_num_modules - length(module_to_be_deleted);
		fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d ',total_num_modules,sum(isactive_module==1));
		iter = iter + 1; % iter is only serving for the drawnow in the following lines
		if round(iter/10)==iter/10
			drawnow;
		end
	end
	fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%7d %7d',total_num_modules,sum(isactive_module==1));drawnow;
	fprintf('\n')
	idx_new = standardize_idx(idx);
	return    

%% calc_module_median:
function [module_median] = calc_module_median(data, idx_samples)

	possible_clusters = setdiff(unique(idx_samples),0);
	module_median = zeros(size(data,1), length(possible_clusters));
	for ii=1:length(possible_clusters)
		ind = find(idx_samples==possible_clusters(ii));
		module_median(:, ii) = median(data(:,ind),2);
	end

	return