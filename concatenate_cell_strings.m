function [outStr] = concatenate_cell_strings(repStr, intArray)
% concatenate_cell_strings.m
% Author: 		Will Macnair 
% Date: 		15/11/13
% Purpose:		
% Relates to:	

intStr = int2str(intArray);
[nn ll] = size(intStr);
labs = repmat({repStr},nn,1); 
ints = mat2cell(intStr, repmat(1,nn,1), ll);
ints = cellfun(@strtrim, ints, 'Unif', false);
allStr = horzcat(labs, ints);
outStr = arrayfun(@(r) [allStr{r,:}], 1:size(allStr,1),'Unif',false);

end
