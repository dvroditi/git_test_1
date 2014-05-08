function [] = check_treeSNE_parameters(treeSNE_parameters)

	if ~isempty(intersect(treeSNE_parameters.used_markers, treeSNE_parameters.extra_markers))
		error('used_markers and extra_markers are not disjoint');
	end