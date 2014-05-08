function for_3_Prepare_DataForSPADE_treeSNE(treeSNE_parameters)

load(treeSNE_parameters.PooledDownsampledDataFilename);

% downsamples data down to manageable value
if size(data,2) <= treeSNE_parameters.num_cells
    keep_ind = 1:size(data,2);
else
    keep_ind = sort(randsample(1:size(data,2), treeSNE_parameters.num_cells));
end

data = data(:,keep_ind);

downsample_factor = numel(keep_ind)/sum(keep_ind);

local_density = local_density(keep_ind)*downsample_factor;
% local_density = local_density(keep_ind);

save(treeSNE_parameters.FlowSPADE_input_filename,'data', 'local_density', 'marker_names', 'used_markers', 'extra_markers');

display(' ')
display(['Prepared SPD input file: ',treeSNE_parameters.FlowSPADE_input_filename,' (', num2str(size(data,2)),' cells)'])
display(' ')
