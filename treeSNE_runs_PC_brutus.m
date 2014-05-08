%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs to establish a robust MST distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prostate

% run run

% define run: tiny run to check things work; takes <1 minute to run
cluster = parcluster('BrutusLSF8h')
matlabpool(cluster,36)

% add appropriate paths
% addpath(genpath('~/Code/'));
% addpath(genpath('~/Packages/'));

% matlabpool close force local; % If you want to parallel toolbox on cluster
% matlabpool; % matlab parallel computing toolbox

treeSNE_parameters.saveStem					= 'treeSNE_test_prostate_phen_2';

% outDir										='/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/'; % Directory to output results
% 
% outDir 										= '~/Outputs/PC/'

outDir 										= '~/tsne/Code/Outputs/PC/'
% Checks if output folder already exists, otherwise creates it

folder_path=[outDir treeSNE_parameters.saveStem];

does_folder_exist 							= exist(folder_path);

% If folder does not exist create directory, otherwise abort and ask for a new folder name
	if does_folder_exist == 0 
	mkdir(outDir,treeSNE_parameters.saveStem);
	else
	error('folder already exists, rename saveStem')
	end


treeSNE_parameters.outputDir				= [outDir treeSNE_parameters.saveStem];

treeSNE_parameters.dataDir					= '~/data/';

% treeSNE_parameters.outputDir				= ['~/Outputs/08_treeSNE/Bodenmiller_PC/' treeSNE_parameters.saveStem];

% treeSNE_parameters.outputDir				= ['/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/'];

treeSNE_parameters.used_markers				= [{'CD10'}, {'CD107a'}, {'CD193'}, {'CD20'}, {'CD3'}, {'CD33'}, {'CD44'}, {'CD45'}, {'CD61'}, {'CD68'}, {'CD71'}, {'CD82'}, {'CD90'}, {'E-cadherin'}]';


treeSNE_parameters.extra_markers			= [{'CAH9'}, {'CK6'}, {'Cleaved Cas3'}, {'CXCR4'}, {'ER'}, {'FAP'}, {'Her2'}, {'KI-67'}, {'pAkt'}, {'pErk1/2'}, {'pNFkB'}, {'pPlcg2'}, {'pRb'}, {'pS6'}, {'pSTAT5'}, {'VEGFR2'}]'; % In case we want to only cluster on a few markers we can add the remaining markers here. If we cluster on all markers, then still add an empty "extra markers" vector here.

% treeSNE_parameters.filenames				= [{'PCG1490-1.fcs'},{'PCG1490-2.fcs'},{'PCG1490-3.fcs'},{'PCG1490-4.fcs'},{'PCG1496-1.fcs'},{'PCG1496-2.fcs'},{'PCG1496-3.fcs'},{'PCG1496-4.fcs'}]'; %[{},...,{}]'

treeSNE_parameters.filenames				= [{'PCG1490-1_fixed.fcs'}, {'PCG1490-2_fixed.fcs'}, {'PCG1490-3_fixed.fcs'}, {'PCG1490-4_fixed.fcs'}, ...
											{'PCG1496-1_fixed.fcs'}, {'PCG1496-2_fixed.fcs'}, {'PCG1496-3_fixed.fcs'}, {'PCG1496-4_fixed.fcs'}]';

treeSNE_parameters.file_annot				= [{'PC_prelim'}]'; % renaming file names for graphs

treeSNE_parameters.fcs_cell_limit			= 10000; % How many cell events from each sample (same as file)

treeSNE_parameters.target_prctile			= 3500; % PERCENTILE of cells to be used if value is < 100; NUMBER of cells if value > 100; done for each file

treeSNE_parameters.exclude_prctile			= 1; % Percentile of cells to be excluded for initial downsampling (thresholding); this is done for each file

treeSNE_parameters.num_cells				= 10000; % Total # of cells actually used from pooled samples (all files); uniformally sampled; SPADE says it breaks down if number is larger than 50K

treeSNE_parameters.num_target_clusters		= 2000; % Self-explanatory

treeSNE_parameters.nDims					= 2; % Dimensions for T-SNE projections

treeSNE_parameters.FlowSPADE_input_filename	= 'DataForSPADE.mat'; % file to hold pooled data

treeSNE_parameters.pool_flag				= true; % To turn on parallel toolbox or not (true or false)

treeSNE_parameters.refSeed					= 12345; % Random number generator seed for sampling reference cells

treeSNE_parameters.num_sample_cells			= treeSNE_parameters.num_target_clusters; % # reference cells, set to be the same as the number of clusters

treeSNE_parameters.nRuns					= 2; % # cluster runs for the same initial set of reference cells; to be averaged later

treeSNE_parameters.arcsinh_cofactor			= 5; % mass cytometry pre-processing (squishes large values since marker readout results in very large values); >5 is treshold for large

treeSNE_parameters.is_normalize				= 0; % SPADE normalization across markers or samples, not sure (?); Turned off for the moment

treeSNE_parameters.kernel_width_para		= 5; % Defines neighborhood size for local density estimate 

% run run

treeSNE(treeSNE_parameters);

treeSNE_robust_distance_estimator_define_ref_cells(treeSNE_parameters);

parfor i = 1:treeSNE_parameters.nRuns
	treeSNE_robust_distance_estimator(treeSNE_parameters, i);
end % parfor

treeSNE_robust_distance_aggregator(treeSNE_parameters);


