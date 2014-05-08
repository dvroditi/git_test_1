%% Initialize/Load data
addpath(genpath('lib'));

clear all;

saveStem='prostate_run_split_large_mean'; %original saveStem from treeSNE_runs plus include _mean
outputDir = '/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/prostate_run_split_large/'; %directory where treeSNE outputs can be found
A = importdata(fullfile(outputDir, [saveStem '_clusters.txt']))


% order rows according to cluster assignement 
[A_sorted,A_ids]=sortrows(A.data);

[unique_clusterID,row_clusterID]=unique(A_sorted(:,1));

% load('data.mat');

% Data ... Time X species X trajectories
% SpeciesCount ... Number of species
% SpeciesDistance ... species x species matrix, true distance
% SpeciesDistanceLabels ... species x species cell array of strings, true distance
% SpeciesNames ... Names of species
% TimepointCount ... Number of time points
% Timepoints ... In seconds
% TrajectoryCount ... Number of Trajectories

% TimepointIdxAnalyze = 41:61;
% TimepointIdxAnalyzeCount = length(TimepointIdxAnalyze);

Lambda = 1e6;
markers=30;
% sigmaHat=;%data matrix
SpeciesNames=[{'ER'}, {'pSTAT5'}, {'VEGFR2'}, {'CD90'}, {'CD71'}, {'CD68'}, {'CD3'}, {'CD33'}, {'CD20'}, {'pNFkB'}, {'Her2'}, ...
{'CD107a'}, {'pAkt'}, {'pErk1/2'}, {'CD10'}, {'CD82'}, {'CD44'}, {'CXCR4'}, {'CK6'}, {'CAH9'}, {'E-cadherin'}, ...
{'KI-67'}, {'pPlcg2'}, {'FAP'}, {'pS6'}, {'CD61'}, {'Cleaved Cas3'}, {'CD193'}, {'pRb'}, {'CD45'}]; %cell array with marker values
% SpeciesDistanceLabels=;
%% Infer GMRF
tic
% sigmaInvs = zeros(SpeciesCount, SpeciesCount, TimepointIdxAnalyzeCount);

sigmaInvs = zeros(markers, markers,length(row_clusterID));

for i=1:length(row_clusterID)-1
	% marker_avg(i,:)=mean(A_sorted(row_clusterID(i):row_clusterID(i+1)-1,:));

% for tIdx = 1:TimepointIdxAnalyzeCount
    sigmaHat = cov(A_sorted(row_clusterID(i):row_clusterID(i+1)-1,:));
	sigmaInvs(:,:,i) = covsel(sigmaHat,Lambda,1e-6);
% end
end
marker_avg(length(row_clusterID),:)=mean(A_sorted(row_clusterID(end):end,:));

toc

% Save heatmap for each timepoint to 'plots/FILENAME_<timepointNo>.png'
% heatmapNormalizedMatrices(sigmaInvs, SpeciesNames, SpeciesDistanceLabels, 'Admm', 'plots/admm2_', Timepoints(TimepointIdxAnalyze));
% heatmapNormalizedMatrices(sigmaInvs, SpeciesNames, SpeciesDistanceLabels, 'Admm', 'plots/admm2_');