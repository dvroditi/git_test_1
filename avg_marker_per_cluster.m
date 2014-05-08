% This function will sort through reference cells belonging to different clusters 
% and calculate the average marker values in each cluster

saveStem='prostate_run_split_large_mean'; %original saveStem from treeSNE_runs plus include _mean
outputDir = '/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/prostate_run_split_large/'; %directory where treeSNE outputs can be found
A = importdata(fullfile(outputDir, [saveStem '_clusters.txt']))


% order rows according to cluster assignement 
[A_sorted,A_ids]=sortrows(A.data);

% Determine the row number where cluster id changes:

[unique_clusterID,row_clusterID]=unique(A_sorted(:,1));

for i=1:length(row_clusterID)-1
	marker_avg(i,:)=mean(A_sorted(row_clusterID(i):row_clusterID(i+1)-1,:));
end
marker_avg(length(row_clusterID),:)=mean(A_sorted(row_clusterID(end):end,:));

RowLabelsValue=[{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'}];

ColumnLabelsValue=[{'ER'}, {'pSTAT5'}, {'VEGFR2'}, {'CD90'}, {'CD71'}, {'CD68'}, {'CD3'}, {'CD33'}, {'CD20'}, {'pNFkB'}, {'Her2'}, ...
{'CD107a'}, {'pAkt'}, {'pErk1/2'}, {'CD10'}, {'CD82'}, {'CD44'}, {'CXCR4'}, {'CK6'}, {'CAH9'}, {'E-cadherin'}, ...
{'KI-67'}, {'pPlcg2'}, {'FAP'}, {'pS6'}, {'CD61'}, {'Cleaved Cas3'}, {'CD193'}, {'pRb'}, {'CD45'}]';

marker_avg=marker_avg(:,2:end);

for ii=1:size(marker_avg,1)
  marker_avg(ii,:)=(marker_avg(ii,:)-min(marker_avg,[],1))./max(marker_avg,[],1);
end 

% marker_avg=(marker_avg-min(marker_avg))./max(marker_avg);

% figure(1)

hmo=HeatMap(marker_avg,'RowLabels', RowLabelsValue,'ColumnLabels', ColumnLabelsValue, 'Colormap',redbluecmap);

% plot(hmo);

addTitle(hmo, 'Average marker value per cluster');


saveas(gcf,'marker_HeatMap.eps', 'psc2');
% print -dbmp marker_HeatMap.bmp
