%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs to establish a robust MST distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define run: tiny run to check things work; takes <1 minute to run
matlabpool;
treeSNE_parameters.saveStem					= 'treeSNE_test';
treeSNE_parameters.dataDir					= '\\d\dfs\Groups\biol\sysbc\claassen\macnairw\Data\Cytobank\experiment_4955_files';
treeSNE_parameters.outputDir				= ['\\d\dfs\Groups\biol\sysbc\claassen\macnairw\Data\Cytobank\' treeSNE_parameters.saveStem];
treeSNE_parameters.used_markers				= [{'Vimentin'}, {'N-Cadherin'}, {'CA125'}, {'CD90'}, {'CD133'}, {'CD10'}, ...
											{'E-Cadherin'}, {'Endoglin'}, {'CD24'}, {'CD44'}, {'CD13'}, {'Mesothelin'}]';
treeSNE_parameters.extra_markers			= [{'pAMPK'}, {'pATM'}, {'pH2AX'}, {'Cyclin-B1'}, {'pNFkB'}, {'pBCL2'}, {'pERK'}, ...
											{'Ki67'}, {'pCHK2'}, {'pSTAT3'}, {'SNAIL'}, {'pAKT'}, {'SOX2'}, {'c-MYC'}, ...
											{'pSTAT5'}, {'pRB'}, {'PAX8'}, {'pCHK1'}, {'NonP-b-cat'}, {'pS6'}, ...
											{'pCREB'}, {'Tot-p53'}, {'pHH3'}]';
treeSNE_parameters.filenames				= [{'20140108_OC2_Wicking_A1-OC-3_PT.fcs'}, {'20140108_OC2_Wicking_A2-U937.fcs'}, ...
											{'20140108_OC2_Wicking_A3-KURAMOCHI.fcs'}, {'20140108_OC2_Wicking_A4-PBMC1.fcs'}]';
treeSNE_parameters.file_annot				= [{'OC-3_PT'}, {'U937'}, {'KURAMOCHI'}, {'PBMC1'}]';
treeSNE_parameters.fcs_cell_limit			= 1000;
treeSNE_parameters.exclude_prctile			= 1;
treeSNE_parameters.target_prctile			= 5000;
treeSNE_parameters.num_cells				= 1000;
treeSNE_parameters.num_target_clusters		= 300;
treeSNE_parameters.nDims					= 2;
treeSNE_parameters.FlowSPADE_input_filename	= 'DataForSPADE.mat';
treeSNE_parameters.pool_flag				= true;
% treeSNE_parameters.pool_flag				= false;
treeSNE_parameters.refSeed					= 123;
treeSNE_parameters.nRuns					= 3;
treeSNE_parameters.num_sample_cells			= treeSNE_parameters.num_target_clusters;
treeSNE_parameters.arcsinh_cofactor			= 5;
treeSNE_parameters.is_normalize				= 0;
treeSNE_parameters.kernel_width_para		= 5;

% run run
treeSNE(treeSNE_parameters);
treeSNE_robust_distance_estimator_define_ref_cells(treeSNE_parameters);

parfor i = 1:treeSNE_parameters.nRuns
	treeSNE_robust_distance_estimator(treeSNE_parameters, i);
end % parfor

treeSNE_robust_distance_aggregator(treeSNE_parameters);


% define run: larger run
matlabpool;
treeSNE_parameters.saveStem					= 'treeSNE_robust_distance_phenotype_only';
treeSNE_parameters.dataDir					= '\\d\dfs\Groups\biol\sysbc\claassen\macnairw\Data\Cytobank\experiment_4955_files';
treeSNE_parameters.outputDir				= ['\\d\dfs\Groups\biol\sysbc\claassen\macnairw\Data\Cytobank\' treeSNE_parameters.saveStem];
treeSNE_parameters.used_markers				= [{'Vimentin'}, {'N-Cadherin'}, {'CA125'}, {'CD90'}, {'CD133'}, {'CD10'}, ...
											{'E-Cadherin'}, {'Endoglin'}, {'CD24'}, {'CD44'}, {'CD13'}, {'Mesothelin'}]';
treeSNE_parameters.extra_markers			= [{'pAMPK'}, {'pATM'}, {'pH2AX'}, {'Cyclin-B1'}, {'pNFkB'}, {'pBCL2'}, {'pERK'}, ...
											{'Ki67'}, {'pCHK2'}, {'pSTAT3'}, {'SNAIL'}, {'pAKT'}, {'SOX2'}, {'c-MYC'}, ...
											{'pSTAT5'}, {'pRB'}, {'PAX8'}, {'pCHK1'}, {'NonP-b-cat'}, {'pS6'}, ...
											{'pCREB'}, {'Tot-p53'}, {'pHH3'}]';
treeSNE_parameters.filenames				= [{'20140108_OC2_Wicking_A1-OC-3_PT.fcs'}, {'20140108_OC2_Wicking_A2-U937.fcs'}, ...
											{'20140108_OC2_Wicking_A3-KURAMOCHI.fcs'}, {'20140108_OC2_Wicking_A4-PBMC1.fcs'}, ...
											{'20140108_OC2_Wicking_A5-ES2.fcs'}, {'20140108_OC2_Wicking_B1-OVCAR4.fcs'}, ...
											{'20140108_OC2_Wicking_B2-COV362.fcs'}, {'20140108_OC2_Wicking_B3-OC-5_OM.fcs'}, ...
											{'20140108_OC2_Wicking_B4-COAV3.fcs'}, {'20140108_OC2_Wicking_B5-OC4_PT.fcs'}, ...
											{'20140108_OC2_Wicking_C1-OC4_OM.fcs'}, {'20140108_OC2_Wicking_C2-OVSAHO.fcs'}, ...
											{'20140108_OC2_Wicking_C3-TYK-NU_Pt_res.fcs'}, {'20140108_OC2_Wicking_C4-OC5_PT.fcs'}, ...
											{'20140108_OC2_Wicking_C5-OVCAR3.fcs'}, {'20140108_OC2_Wicking_D1-TYK-NU.fcs'}, ...
											{'20140108_OC2_Wicking_D2-PBMC2.fcs'}, {'20140108_OC2_Wicking_D3-OC-3_OM.fcs'}, ...
											{'20140108_OC2_Wicking_D4-MCF7.fcs'}, {'20140108_OC2_Wicking_D5-OVKATE.fcs'}]';
treeSNE_parameters.file_annot				= [{'OC-3_PT'}, {'U937'}, {'KURAMOCHI'}, {'PBMC1'}, {'ES2'}, {'OVCAR4'}, {'COV362'}, {'OC-5_OM'}, {'COAV3'}, ...
											{'OC4_PT'}, {'OC4_OM'}, {'OVSAHO'}, {'TYK-NU_Pt_res'}, {'OC5_PT'}, {'OVCAR3'}, {'TYK-NU'}, {'PBMC2'}, {'OC-3_OM'}, ...
											{'MCF7'}, {'OVKATE'}]';
% treeSNE_parameters.filenames				= [{'20140108_OC2_Wicking_A1-OC-3_PT.fcs'}, {'20140108_OC2_Wicking_A2-U937.fcs'}, ...
% 											{'20140108_OC2_Wicking_A3-KURAMOCHI.fcs'}, {'20140108_OC2_Wicking_A4-PBMC1.fcs'}]';
% treeSNE_parameters.file_annot				= [{'OC-3_PT'}, {'U937'}, {'KURAMOCHI'}, {'PBMC1'}]';
treeSNE_parameters.fcs_cell_limit			= 5000;
treeSNE_parameters.exclude_prctile			= 1;
treeSNE_parameters.target_prctile			= 20000;
treeSNE_parameters.num_cells				= 10000;
treeSNE_parameters.num_target_clusters		= 1000;
treeSNE_parameters.nDims					= 2;
treeSNE_parameters.FlowSPADE_input_filename	= 'DataForSPADE.mat';
treeSNE_parameters.pool_flag				= true;
% treeSNE_parameters.pool_flag				= false;
treeSNE_parameters.refSeed					= 123;
treeSNE_parameters.nRuns					= 100;
treeSNE_parameters.num_sample_cells			= treeSNE_parameters.num_target_clusters;
treeSNE_parameters.arcsinh_cofactor			= 5;
treeSNE_parameters.is_normalize				= 0;
treeSNE_parameters.kernel_width_para		= 5;

% run run
treeSNE(treeSNE_parameters);
treeSNE_robust_distance_estimator_define_ref_cells(treeSNE_parameters);

parfor i = 1:treeSNE_parameters.nRuns
	treeSNE_robust_distance_estimator(treeSNE_parameters, i);
end % parfor

treeSNE_robust_distance_aggregator(treeSNE_parameters);
