% Cross validation for finding best lambda

%% Initialize/Load data
addpath(genpath('lib'));
clear;
load('data.mat');
% TrajectoryCount = 5;
% SpeciesCount = 6;
% spidx = [1 3 4 8 10 14];
% SpeciesNames = SpeciesNames(spidx); 
% DataQT(:,:,TrajectoryCount+1:end)=[];
% DataQT(:,setdiff(1:14, spidx),:)=[];
% Data(:,:,TrajectoryCount+1:end)=[];
% Data(:,setdiff(1:14, spidx) ,:)=[];
% SpeciesDistanceLabels(setdiff(1:14, spidx),:) = [];
% SpeciesDistanceLabels(:,setdiff(1:14, spidx)) = [];

% Data ... Time X species X trajectories
% SpeciesCount ... Number of species
% SpeciesDistance ... species x species matrix, true distance
% SpeciesDistanceLabels ... species x species cell array of strings, true distance
% SpeciesNames ... Names of species
% TimepointCount ... Number of time points
% Timepoints ... In seconds
% TrajectoryCount ... Number of Trajectories

% Parameter for computation
% TimepointIdxAnalyze = 1:length(Timepoints);
% TimepointIdxAnalyzeCount = length(TimepointIdxAnalyze);
Lambdas = [0 logspace(-9,1,10)];
CVFold = 5;
relTol = 1e-4;
TrajectoryCount=30;

%% Infer GMRF

% Prepare training and test sets
if (mod(TrajectoryCount, CVFold))
    error('TrajectoryCount is not a multiple of CrossValidationFold.');
end
CVTestSets = reshape(randperm(TrajectoryCount), [TrajectoryCount/CVFold CVFold ]);
CVTrainingSets = zeros(TrajectoryCount-size(CVTestSets,1), CVFold);
for cvIdx = 1:CVFold
    CVTrainingSets(:,cvIdx) = setdiff(1:TrajectoryCount,CVTestSets(:,cvIdx));
end

% Store best lambda for each time point
% SelectedLambdas = zeros(TimepointIdxAnalyzeCount, 1);

% sigmaInvs = zeros(SpeciesCount,SpeciesCount,TimepointIdxAnalyzeCount);

sigmaInvs = zeros(markers,markers);
tic
% for tIdx = 1:TimepointIdxAnalyzeCount
     tpData = squeeze(DataQT);
%     tpData = squeeze(Data(TimepointIdxAnalyze(tIdx),:,:));
    
    % Intermediate results for lambda
    logLikLambdas = zeros(length(Lambdas),1);
    for lIdx = 1:length(Lambdas)
        lambda = Lambdas(lIdx); % lambda for current iteration
        
        % LogLik for each cross validation
        logLikCVs = zeros(CVFold, 1);
        for cvIdx = 1:CVFold
            trainingData = tpData(:,CVTrainingSets(:,cvIdx));
            testData = tpData(:,CVTestSets(:,cvIdx));
            
            % Compute empirical covariance for trainingset
            sigmaHat = cov(trainingData');
            
            % Compute empirical covariance for testset
            sigmaHatTest = cov(testData');
            
            % Run optimization with training data
            sigmaInv = covsel(sigmaHat, lambda, relTol);

            % Evalulate logLik on test set
            logLikCVs(cvIdx) = log(det(sigmaInv)) - trace(sigmaHatTest*sigmaInv);
        end
        logLikLambdas(lIdx) = mean(logLikCVs);
    end
    % Selected best lambda
    [bestLL, bestLLIdx] = max(logLikLambdas); 
    % SelectedLambdas(tIdx) = Lambdas(bestLLIdx);
    
    % Recompute/store sigma^-1 with best lambda and entire dataset
    sigmaHat = cov(tpData');
    sigmaInvs = covsel(sigmaHat, Lambdas(bestLLIdx), relTol);
% end
toc

% Plot time course of best lambdas
bar(log10(SelectedLambdas)); title('Lambdas selected by cross validation'); xlabel('time point index'); ylabel('log10(lambda)');

% Generate pdf with GMRFs
heatmapNormalizedMatrices(sigmaInvs, SpeciesNames, SpeciesDistanceLabels, 'ADMM time=', 'plots/admmcv1000_', Timepoints(TimepointIdxAnalyze));
