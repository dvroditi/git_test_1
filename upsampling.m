% go to directory with treeSNE files

% function []=upsampling(a,b)

% (pooled_data,mean_clusters)

% addpath('/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/treeSNE_test_prostate_phen_pool_mex');

addpath('/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/treeSNE_run_Will_server');

saveStem = 'treeSNE_run_Will_server'; %for consistency, use the same as treeSNE_runs_... file.

outDir ='/Users/delaura/Documents/MATLAB/Code/treeSNE/PC/'; % Directory to outpout results


load PooledData.mat; % Uploads a data matrix of all cells pooled from all 
                     % files used for treeSNE and all marker values collected
                     % Also uploads a one-dimensional array of indexes of 
                     % which markers were used for clustering - used_markers 
                     % as well as an array of extra markers and a cell array
                     % of marker_names. 
pooled_cells=data;

pooled_data_used_markers=pooled_cells(used_markers, :); % an array containing marker values for 
                                                        % all pooled data for only markers used 
                                                        % for clustering

mean_clusters=importdata('prostate_run_split_large_mean_clusters.txt'); % Struct containing data: [double]; 
                                                                                   % textdata: {cell}; colheaders: {cell} 

mean_cluster_clusterID_used_and_extra_markers=mean_clusters.data; % data is the mean cluster value (similar to a 
                                                                 % reference cell) and corresponding marker values
                                                                 % as well as MASTER cluster assingment

cl_hd=mean_clusters.colheaders; % cell array containing "used markers"
                                %"extra markers" and "cluster" 

used_marker_labels=marker_names(used_markers);

[A,B]=ismember(cl_hd,used_marker_labels); % A is an array of logicals so we  know the index location of the 
                                          % markers that match between the column header of the mean_clusters 
                                          % file and the used markers from PooledData

A(12)=1;A(15)=1;A(23)=1;A(22)=1;A(28)=1; % Be careful! If names are not exactly the same, you will not get 
         % a match, this is particularly important in the case of E-cadherin 
         % (versus E.cadherin) which seems to be spelled differently in          
         % different places. Here I just force E-cadherin to be a match since 
         % there was a difference in spelling from one file to the next
           
A(1)=1; %This is to include the first element of column header (cl_hd) which is the cluster ID which isnt' in u

mean_cluster_used_markers=mean_cluster_clusterID_used_and_extra_markers(:,A); % marker values of mean clusters ONLY 
                                                         % for markers used in clustering 
                                                         
mean_cluster_used_markers_noID=mean_cluster_used_markers(:,(2:end)); % Before calculating the distance of cells to cluster 
                                                                     % means we need to exclude the first column containing 
                                                                     % "Master" cluster ID
                                                                                                          
mean_cluster_used_markers_noID=mean_cluster_used_markers_noID'; % for "distance" function, markers need to be down the rows
mean_cluster_used_markers=mean_cluster_used_markers';

E=distance(pooled_data_used_markers,mean_cluster_used_markers_noID);                                                                    

[minColVal, minColIdx] = min(E,[],2); % min is taken across columns since we want the min values across mean clusters for each individual cell

cluster_ids=mean_cluster_used_markers(1,:);
mean_cl=1:size(mean_cluster_clusterID_used_and_extra_markers,1);
mean_cl=mean_cl';
clusterID_meanCluster=[cluster_ids' mean_cl];
data_mean_cluster_location=[data; minColIdx']; % The last row of this matrix (minColIdx') contains the mean cluster closest to a particular corresponding cell on the column

% Now we adjust the minColIdx array to represent the "MASTER" cluster to
% which each individual cell belongs to
num_clusters=max(cluster_ids);
% mean_cl_ids=zeros(num_clusters,) % second argument has varied size....
minColIdx_temp=minColIdx;
for i=1:num_clusters
    
    mean_cl_ids=find(cluster_ids==i);
    index_temp= find(ismember(minColIdx,mean_cl_ids)); % minColIdx and minColIdx_temp are initially the same.
                                                       % We look for the indexes on minColIdx and modify minColIdx_temp
                                                       % so mean cluster values that have already been converted to a 
                                                       % Master cluster don't get reassigned to another Master cluster. 
                                                       % Ex: mean cluster 2 belongs to Master cluster 3, and therefore 
                                                       % will be reassigned the value 3, but mean cluster 
                                                       % 3 belongs to Master cluster 5. this can cause the misassigment of
                                                       % mean cluster 2 to Master cluster 5 therefore we check
                                                       % for indexes in the unmodified file minColIdx
 
    minColIdx_temp(index_temp)=i;
    
end

data_Master_cluster_location=[data_mean_cluster_location; minColIdx_temp']; %data file containing individual cells and their Master cluster assignment. 

for i=1:num_clusters
    
a=find(minColIdx_temp==i); % find indexes for cells in cluster i

data_in_upsampled_cluster=data(:,a);

saveFilename = [saveStem '_indexes_upsampled_cluster_' num2str(i) '.txt'];
  
saveFilepath=fullfile(outDir,saveStem,saveFilename);

save (saveFilepath, 'a', '-ascii', '-tabs');

end

% set up module save
	
%     txtfile=['File' num2str(i) '.txt']; % Need to convert index to string
% 	header = extra_marker_names';
% 	saveData = mean_extra_location;
% 	% save density
% 	save_txt_file(saveFilename, header, saveData);

% function [] = save_txt_file(saveFilename, header, saveData)
% 	if size(header,2) ~= size(saveData,2)
% 		error(['Problem saving ' saveFilename ': header and data are not compatible lengths.']);
% 	else
% 		disp(['Saving file ' saveFilename]);
% 	end
% 
% 	headerSpec = ['%s' repmat('\t%s', 1, size(header,2) -1) '\n'];
% 	dataSpec = ['%4.4f' repmat('\t%4.4f', 1, size(saveData,2) -1) '\n'];;
% 	fid = fopen(saveFilename, 'w');
% 	fprintf(fid, headerSpec, header{:});
% 	for ii = 1:size(saveData,1)
% 		fprintf(fid, dataSpec, saveData(ii,:));
% end
% fclose(fid);