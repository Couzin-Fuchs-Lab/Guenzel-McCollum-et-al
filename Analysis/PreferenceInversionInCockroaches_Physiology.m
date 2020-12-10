%% PreferenceInversionInCockroaches_Physiology
% 
% Summary:
% In social species, decision-making is both influenced by, and in turn
% influences, the social context. This reciprocal feedback introduces 
% coupling across scales, from the neural basis of sensing, to individual
% and collective decision-making. Here, we adopt an integrative approach 
% investigating decision-making in dynamical social contexts. When choosing
% shelters, isolated cockroaches prefer vanillin-scented (food-associated) 
% shelters over unscented ones, yet in groups, this preference is inverted.
% We demonstrate that this inversion can be replicated by replacing the 
% full social context with social odours: presented alone food and social 
% odours are attractive, yet when presented as a mixture they are avoided. 
% Via antennal lobe calcium imaging, we show that neural activity in 
% vanillin-responsive regions reduces as social odour concentration 
% increases. Thus, we suggest that the mixture is evaluated as a distinct 
% olfactory object with opposite valence, providing a mechanism that 
% would naturally result in individuals avoiding what they perceive as 
% recently exploited resources.
%
% This script represents the pipeline used to analyse calcium-imaging data.
% It uses baseline-subtracted Delta R340/380 2D fluorescence time courses
% and applies a k-medoids cluster analysis on active pixels 
% (fluorescence > 0.01 for at least one of the five odours).
%
% (MATLAB 2020a)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

% Add path to save figures
addpath(genpath('altmany-export_fig'))

% Create folder to save plots
if ~isdir('FIG\raw_imaging'); mkdir('FIG\raw_imaging'); end


%% Settings

% Get where the data is stored. Note, this is not the raw data as one gets
% it from the experiments (pst files), but slightly processed secondary
% data. What has been done is performing the ratio of 340/380 nm videos and
% subtraction of the baseline.
SET.Path2Data = 'DATA\Physiology\';
SET.FileNames = {...
    '181217_mp03_dff.mat';...
    '181219_mp04_dff.mat';...
    '190430_ip23_dff.mat';...
    '190508_ip28_dff.mat';...
    '190523_ip33_dff.mat';...
    '191119_ip40_dff.mat';...
    '191120_ip43_dff.mat';...
    '191121_mp43_dff.mat';...
    '191121_mp44_dff.mat';...
    '191122_mp46_dff.mat';...
    };

% Set smothing window
SET.Smooth = [4 4];

% Set whether to use responding pixels only
SET.RespondingOnly = 1;

% Set threshold for active pixel
SET.threshold = 0.01;

% Set parameters for PCA-based feature reduction
SET.PCA_VarExplained = 95; %use PCs to get the variance explained

% Set kmedoids parameter
SET.kmedoids.nCluster = 25; % 5stimuli ^ 2ways of response = 25 possibilities
SET.kmedoids.nIter = 10000;
SET.kmedoids.nRep = 100;
SET.kmedoids.OnlinePhase = 'on';
SET.kmedoids.UseParallel = true;
SET.kmedoids.Dist = 'sqeuclidean';
SET.kmedoids.Start = 'cluster';
SET.kmedoids.EvalCriterion = 'silhouette';

% Set odor names we will use
SET.OdorNames = {...
    'van2';...
    'fec3van2';...
    'fec2van2';...
    'fec2van3';...
    'fec2'};

% Color code for different conditions
SET.Color = [...
    112   174   224;... Van2
    172   215   232;... Van2Fe3
    188   205   188;... Van2Fe2
    232   215   172;... Van3Fe2
    224   174   112;... Fe2
    ]/255;

%% Load, restructure, smooth data, followed by subtracting the response to mineral oil and averaging repetitions

% Iterate over all files
for iFile = 1:length(SET.FileNames)
    
    
    % ---(1)---------------------------------------------------------------
    % Load file
    currData = load([SET.Path2Data, SET.FileNames{iFile}]);
    
    
    % ---(2)---------------------------------------------------------------
    % Iterate over all tested stimuli, but only keep those we are
    % interested in
    for iOdor = 1:length(currData.odor_names)
        % Adjust odor names
        switch currData.odor_names{iOdor}
            case 'fec-2'
                currData.odor_names{iOdor} = 'fec2';
            case 'van-2'
                currData.odor_names{iOdor} = 'van2';
            case 'fec-3+van-3'
                currData.odor_names{iOdor} = 'fec3van3';
            case 'fec-2+van-3'
                currData.odor_names{iOdor} = 'fec2van3';
            case 'fec-3+van-2'
                currData.odor_names{iOdor} = 'fec3van2';
            case 'fec-2+van-2'
                currData.odor_names{iOdor} = 'fec2van2';
        end%switch
        % Only take interesting stimuli: van2, fec3van2, fec2van2, fec2van3, fec2
        if sum(strcmp(SET.OdorNames, currData.odor_names{iOdor}))>0 || strcmp('mol', currData.odor_names{iOdor})
            % --- Iterate over repetitions
            cnt=1;
            for iRep = iOdor:length(currData.odor_names):length(currData.dff)
                currData.Response.(currData.odor_names{iOdor}){cnt} = currData.dff{iRep};
                % Smooth
                currData.SmoothResponse.(currData.odor_names{iOdor}){cnt} = currData.Response.(currData.odor_names{iOdor}){cnt};
                for iFrame = 1:length(currData.vec)
                    currData.SmoothResponse.(currData.odor_names{iOdor}){cnt}(:, :, iFrame) = medfilt2(currData.Response.(currData.odor_names{iOdor}){cnt}(:, :, iFrame), SET.Smooth);
                end
                cnt=cnt+1;
            end%iRep
        end%if interesting stimulus
    end%iOdor
    clear iOdor iRep cnt iFrame
    
    
    % ---(3)---------------------------------------------------------------
    % Subtract response to mol
    % Iterate over repetitions
    for iRep = 1:currData.rep
        % Iterate over odors
        for iOdor = 1:length(SET.OdorNames)
            % Subtract mol-response of the current repetition
            currData.Response.(SET.OdorNames{iOdor}){iRep} = currData.Response.(SET.OdorNames{iOdor}){iRep} - currData.Response.mol{iRep};
            currData.SmoothResponse.(SET.OdorNames{iOdor}){iRep} = currData.SmoothResponse.(SET.OdorNames{iOdor}){iRep} - currData.SmoothResponse.mol{iRep};
        end%iOdor
    end% iRep
    clear iRep iOdor
    
    
    % ---(4)---------------------------------------------------------------
    % Get average response for each odorant
    % Iterate over odors
    for iOdor = 1:length(SET.OdorNames)
        % --- Raw data
        % Preallocation
        temp = zeros(size(currData.Response.mol{1}, 1), size(currData.Response.mol{1}, 2), size(currData.Response.mol{1}, 3));
        % Summ up responses from repetitions
        for iRep = 1:currData.rep
            temp = temp+currData.Response.(SET.OdorNames{iOdor}){iRep};
        end
        clear iRep
        % Divide by number of repetitions
        temp = temp/currData.rep;
        % Save as avgResponse
        currData.avgResponse.(SET.OdorNames{iOdor}) = temp;
        
        % --- Smothed data
        % Preallocation
        temp = zeros(size(currData.SmoothResponse.mol{1}, 1), size(currData.SmoothResponse.mol{1}, 2), size(currData.SmoothResponse.mol{1}, 3));
        % Summ up responses from repetitions
        for iRep = 1:currData.rep
            temp = temp+currData.SmoothResponse.(SET.OdorNames{iOdor}){iRep};
        end
        clear iRep
        % Divide by number of repetitions
        temp = temp/currData.rep;
        % Save as avgResponse
        currData.AvgSmoothResponse.(SET.OdorNames{iOdor}) = temp;
    end%iOdor
    clear iOdor temp
    
    
    % ---(5)---------------------------------------------------------------
    % Combine everything in one data structure
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF = currData.BF_Image;
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).vec = currData.vec;
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse = currData.SmoothResponse;
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).AvgSmoothResponse = currData.AvgSmoothResponse;
    
    % To save memory, stop loop here, take only necessary data, and create
    % the tables used for clustering in a different loop
    clear currData
    
    % Save state
    save('ClusterData', 'DATA', 'SET', '-v7.3')
    
end%iFile
clear iFile

%% Prepare table for clustering
% For this, features are the avg. responses to the different stimuli and
% objects are the responding pixels

% Prespecification
DATA.ClusterResponse_all.mat = [];
DATA.ClusterResponse_all.min_corrcoef = [];

% Waitbar
wb = waitbar(0, 'Please wait...');
% Iterate over all files
for iFile = 1:length(SET.FileNames)
    
    waitbar(iFile/length(SET.FileNames), wb)
    
    % ---(1)---------------------------------------------------------------
    % Get response intervals
    samplenum = str2double(SET.FileNames{iFile}(10:11));
    if samplenum<38
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [27, 29];
    elseif samplenum == 30
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [29, 32];
    elseif samplenum == 32
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [33, 35];
    else
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [32, 34];
    end
    clear samplenum
    
    
    % ---(2)---------------------------------------------------------------
    % Iterate over all pixels and create table with features for each of
    % them
    % Prealloction
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data =           zeros(numel(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF), length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).vec)*length(SET.OdorNames));
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates =    zeros(numel(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF), 2);
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.weakresponse =   zeros(numel(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF), 1);
    % Counter
    cnt = 1;
    % Iterate over all pixels
    for iX = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1)
        for iY = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2)
            % Concatenate respones to the different odors
            temp = [];
            for iOdor = 1:length(SET.OdorNames)
                % Get response
                temp = [temp, squeeze(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).AvgSmoothResponse.(SET.OdorNames{iOdor})(iX, iY, :))'];
                % Check whether the current pixel was responding to, at least,
                % one stimulus
                if nanmean(squeeze(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).AvgSmoothResponse.(SET.OdorNames{iOdor})(iX, iY, SET.Intervals_Stim(1):SET.Intervals_Stim(2)+3))) < SET.threshold
                    WeakFlag(iOdor) = 1;
                else
                    WeakFlag(iOdor) = 0;
                end
            end%iOdor
            DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data(cnt, :) = temp;
            DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(cnt, 1:2) = [iX, iY];
            DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.weakresponse(cnt) = sum(WeakFlag)==length(SET.OdorNames);
            cnt = cnt+1;
            clear  temp WeakFlag iOdor
        end%iY
    end%iX
    clear iX iY cnt
    
    
    % ---(3)---------------------------------------------------------------
    % Kick out not responding pixels
    if SET.RespondingOnly
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.weakresponse==1, :) = [];
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.weakresponse==1, :) = [];
    end
    
    
    % ---(4)---------------------------------------------------------------
    % Apply PCA to reduce feature space
    [~, pcaAnalysis.score, ~, ~, pcaAnalysis.explained, ~] = pca(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data);
    % Get a reduced feature space
    % --- Get the correct number of PCs
    cumulVarExplained = cumsum(pcaAnalysis.explained);
    [~, idx] = min(abs(cumulVarExplained-SET.PCA_VarExplained));
    % --- Get a PC-based table for clustering
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data_PCA = pcaAnalysis.score(:, 1:idx);
    clear pcaAnalysis cumulVarExplained idx
    
    
    % ---(5)---------------------------------------------------------------
    % Now that we have determined the optimal number of clusters, use it to
    % divide the pixels in groups using the kmedoids algorithm
    rng(4243); % For reproducibility
    [DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID, DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids] = ...
        kmedoids(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data_PCA, ...               data
        SET.kmedoids.nCluster, ...                                                                        clusters
        'OnlinePhase', SET.kmedoids.OnlinePhase, ...                                                      Turning on this feature tends to improve the quality of solutions
        'Distance', SET.kmedoids.Dist, ...                                                                'seuclidean' | 'cityblock' | 'minkowski' | 'chebychev' | 'mahalanobis' | 'cosine' | 'correlation' | 'spearman' | 'hamming' | 'jaccard'
        'Options', ...
        statset('MaxIter', SET.kmedoids.nIter, ...                                                        Maximum number of iterations allowed
        'UseParallel', SET.kmedoids.UseParallel), ...                                                     If true, compute in parallel
        'Replicates', SET.kmedoids.nRep, ...                                                              Number of times to repeat clustering using new initial cluster positions.
        'Start', SET.kmedoids.Start...                                                                   'plus' | 'sample' | 'cluster'
        );
    
    
    % ---(6)---------------------------------------------------------------
    % Get each cluster's response to each stimulus
    % Preallocation
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat = zeros(length(SET.OdorNames), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.mol)));
    % Iterate over stimuli
    for iOdor = 1:length(SET.OdorNames)
        % Preallocation
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg = zeros(length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor})), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1));
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std = zeros(length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor})), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1));
        % Iterate over cluster
        for iCluster = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1)
            % Iterate over repetitions
            for iRep = 1:length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor}))
                % Get average response to stimulus
                temp = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor}){iRep}(:, :, SET.Intervals_Stim(1):SET.Intervals_Stim(2)), 3);
                % Get response of current cluster's pixels
                x_coord = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster), 1);
                y_coord = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster), 2);
                px = zeros(length(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster)), 1);
                % Iterate over pixels
                for iPx = 1:length(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster))
                    px(iPx) = temp(x_coord(iPx), y_coord(iPx));
                end%iPx
                % Store for later
                DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).data{iRep, iCluster} = px ;
                DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(iRep, iCluster)  = nanmean(px) ;
                DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std(iRep, iCluster)  = std(px);%/sqrt(length(px)) ;
                DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat(iOdor, iCluster, iRep) = nanmean(px);
                % Tidy up
                clear temp px iPx x_coord y_coord temp
            end%iRep
        end%iCluster
    end%iOdor
    clear iOdor iCluster iRep
    
    
    % ---(7)---------------------------------------------------------------
    % Sort clusters by van2 responsiveness
    [~, idx] = sort(nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.van2.avg), 'descend');
    temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID = zeros(length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID),1);
    for iCluster = 1:length(idx)
        for iOdor = 1:length(SET.OdorNames)
            for iRep = 1:length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor}))
                temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).data{iRep, iCluster} = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).data{iRep, idx(iCluster)};
                temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(iRep, iCluster)  = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(iRep, idx(iCluster));
                temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std(iRep, iCluster)  = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std(iRep, idx(iCluster));
                temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat(iOdor, iCluster, iRep) =                  DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat(iOdor, idx(iCluster), iRep);
            end%iRep
        end%iOdor
        temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == idx(iCluster))) = iCluster;
    end%iCluster
    clear idx iCluster iOdor iRep
    for iOdor = 1:length(SET.OdorNames)
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).data = temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).data;
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg  = temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg;
        DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std  = temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).std;
    end%iOdor
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat = temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat;
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID = temp.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID;
    clear temp iCluster iOdor iRep
    
    
    % ---(8)---------------------------------------------------------------
    % Correlation coefficients
    % Preallocation
    DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.corrcoef = zeros(length(SET.OdorNames), length(SET.OdorNames));
    temp = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat, 3);
    % Iterate over stimuli
    for iOdor_1 = 1:length(SET.OdorNames)
        % Iterate over stimuli
        for iOdor_2 = 1:length(SET.OdorNames)
            % Calculate the correlation coe
            R = corrcoef(temp(iOdor_1, :), temp(iOdor_2, :));
            DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.corrcoef(iOdor_1, iOdor_2) = R(2);
            clear R
        end%iOdor_1
    end%iOdor_2
    clear iOdor_1 iOdor_2 temp
    
    
    % ---(9)---------------------------------------------------------------
    % Combine the responses over all trials
    DATA.ClusterResponse_all.mat = cat(2,DATA.ClusterResponse_all.mat, nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat, 3));
    DATA.ClusterResponse_all.min_corrcoef = [DATA.ClusterResponse_all.min_corrcoef; min(min(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.corrcoef))];
    
    
end%iFile
clear iFile

% ---(10)---------------------------------------------------------------
% Correlation coefficients  over all trials
% Preallocation
DATA.ClusterResponse_all.corrcoef = zeros(length(SET.OdorNames), length(SET.OdorNames));
% Iterate over stimuli
for iOdor_1 = 1:length(SET.OdorNames)
    % Iterate over stimuli
    for iOdor_2 = 1:length(SET.OdorNames)
        % Calculate the correlation coe
        R = corrcoef(DATA.ClusterResponse_all.mat(iOdor_1, :), DATA.ClusterResponse_all.mat(iOdor_2, :));
        DATA.ClusterResponse_all.corrcoef(iOdor_1, iOdor_2) = R(2);
        clear R
    end%iOdor_1
end%iOdor_2
DATA.ClusterResponse_all.min_corrcoef = [DATA.ClusterResponse_all.min_corrcoef; min(min(DATA.ClusterResponse_all.corrcoef))];
clear iOdor_1 iOdor_2

% Save state
save('ClusterData', 'DATA', 'SET', '-v7.3')

% Indicate that clustering is over
close(wb)

%% Create figures for each trial

% Settings aquired post-clustering for examples
SET.Examples.caxis = [-0.0123, 0.0255];
SET.Examples.file_181217_mp03_dff.Cluster = [1, 11, 3, 4, 24, 25];
SET.Examples.file_181217_mp03_dff.ylim = [-0.0114, 0.0288];
SET.Examples.file_181217_mp03_dff.TC_ylim = [-0.015 0.03];
SET.Examples.file_181219_mp04_dff.Cluster = [1, 3, 6, 12, 24, 25];
SET.Examples.file_181219_mp04_dff.ylim = [-0.0215, 0.0315];
SET.Examples.file_181219_mp04_dff.TC_ylim = [-0.02 0.03];

% Iterate over all files
hFig.ClusterCorrelation = figure('Units', 'Normalized', 'Position', [0 0 1 1], 'Color', 'w');
for iFile = 1:length(SET.FileNames)
    
    % ---(1)---------------------------------------------------------------
    % Get response intervals
    samplenum = str2double(SET.FileNames{iFile}(10:11));
    if samplenum<38
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [27, 29];
    elseif samplenum == 30
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [29, 32];
    elseif samplenum == 32
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [33, 35];
    else
        SET.Intervals_BL = [15, 22];
        SET.Intervals_Stim = [32, 34];
    end
    clear samplenum
    
    % ---(Prepare figures)-------------------------------------------------
    hFig.PhysiologyHeatmaps = figure('Units', 'Normalized', 'Position', [0 0 1 1], 'Color', 'w', 'Name', SET.FileNames{iFile}(1:11));
    hFig.ClusterOverview = figure('Units', 'centimeters', 'Position', [0 0 19 20], 'Color', 'w', 'Name', SET.FileNames{iFile}(1:11));
    
    
    % Iterate over all stimuli
    for iOdor = 1:length(SET.OdorNames)
        
        % ---(PhysiologyHeatmaps)------------------------------------------
        figure(hFig.PhysiologyHeatmaps)
        if iOdor == 1 %plot BF
            ax(1) = subplot(1, length(SET.OdorNames)+1, iOdor);
            imagesc(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF')
            axis equal
            xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
            ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
            box off
            set(gca, 'XTick', [], 'YTick', [])
            colormap(ax(1), gray(50))
            cb = colorbar;
            cb.Location = 'southoutside';
            title('BF')
            % Special treatment for examples
            if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
                hold on
                cm = viridis(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1));
                % Iterate over cluster and add boundaries
                for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
                    % Get current cluster's pixels
                    idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
                    % Create a BW image of only the current cluster's
                    % pixels being white
                    BW = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
                    for iPX = 1:length(idx)
                        BW(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 2)) = 1;
                    end%iPX
                    BW = BW';
                    % Draw cluster's outlines
                    [B, L] = bwboundaries(BW, 'noholes');
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'Color', cm(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster), :), 'LineWidth', 1)
                    end%k
                    clear k B L BW idx
                end%iCluster
                clear cm
            end%if special treatment
        end%if BF
        ax(iOdor+1) = subplot(1, length(SET.OdorNames)+1, iOdor+1);
        temp = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).AvgSmoothResponse.(SET.OdorNames{iOdor})(:, :, SET.Intervals_Stim(1):SET.Intervals_Stim(2)), 3);
        imagesc(temp'); clear temp
        axis equal
        xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
        ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
        box off
        set(gca, 'XTick', [], 'YTick', [])
        colormap(ax(iOdor+1), viridis(50))
        % Special treatment for examples
        if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
            caxis(SET.Examples.caxis)
            hold on
            % Iterate over cluster and add boundaries
            for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
                % Get current cluster's pixels
                idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
                % Create a BW image of only the current cluster's
                % pixels being white
                BW = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
                for iPX = 1:length(idx)
                    BW(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 2)) = 1;
                end%iPX
                BW = BW';
                % Draw cluster's outlines
                [B, L] = bwboundaries(BW, 'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                end%k
                clear k B L BW idx
            end%iCluster
            clear cm
        else
            caxis([min(min(DATA.ClusterResponse_all.mat)), max(max(DATA.ClusterResponse_all.mat))])
        end%if
        cb = colorbar;
        cb.Location = 'southoutside';
        title(SET.OdorNames{iOdor})
        
        
        % ---(ClusterOverview: BF)-----------------------------------------
        if iOdor == 1 %plot BF
            figure(hFig.ClusterOverview)
            ax_bf = subplot(7, 7, [iOdor, iOdor+7]);
            imagesc(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF')
            axis equal
            xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
            ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
            box off
            set(gca, 'XTick', [], 'YTick', [])
            colormap(ax_bf, gray(50))
            cb = colorbar;
            cb.Location = 'southoutside';
            title('BF')
            % Special treatment for examples
            if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
                hold on
                cm = viridis(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1));
                % Iterate over cluster and add boundaries
                for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
                    % Get current cluster's pixels
                    idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
                    % Create a BW image of only the current cluster's
                    % pixels being white
                    BW = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
                    for iPX = 1:length(idx)
                        BW(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 2)) = 1;
                    end%iPX
                    BW = BW';
                    % Draw cluster's outlines
                    [B, L] = bwboundaries(BW, 'noholes');
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(boundary(:, 2), boundary(:, 1), 'Color', cm(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster), :), 'LineWidth', 1)
                    end%k
                    clear k B L BW idx
                end%iCluster
                clear cm
            end%if special treatment
            set(gca, 'Fontsize', 3)
        end%if BF
        
        
        % ---(ClusterOverview: clustermap)---------------------------------
        if iOdor == 1 %plot clustermap
            figure(hFig.ClusterOverview)
            ax_cm = subplot(7, 7, [7 14]);
            temp = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
            for iPx = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates, 1)
                temp(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(iPx, 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(iPx, 2)) = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID(iPx);
            end%iPx
            imagesc(temp')
            if min(unique(temp)) == 0
                cm = viridis(length(unique(temp)));
                cm(1, :) = 0;
                colormap(ax_cm, cm)
            else
                colormap(ax_cm, viridis(length(unique(temp))))
            end
            axis equal
            xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
            ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
            box off
            set(gca, 'XTick', [], 'YTick', [])
            cb = colorbar;
            cb.Location = 'southoutside';
            set(cb, 'YTick', linspace(0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1)-0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1)+1), 'YTickLabel', 0:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1))
            title(['Number of active pixel: ', num2str(length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID))])
            set(gca, 'XTick', [], 'YTick', [])
            set(gca, 'Fontsize', 3)
        end%if clustermap
        
        
        % ---(ClusterOverview: heatmaps)-----------------------------------
        figure(hFig.ClusterOverview)
        ax_heat(iOdor) = subplot(7, 7, [iOdor+1, iOdor+1+7]);
        temp = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).AvgSmoothResponse.(SET.OdorNames{iOdor})(:, :, SET.Intervals_Stim(1):SET.Intervals_Stim(2)), 3);
        imagesc(temp'); clear temp
        axis equal
        xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
        ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
        box off
        set(gca, 'XTick', [], 'YTick', [])
        colormap(ax_heat(iOdor), viridis(50))
        % Special treatment for examples
        if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
            caxis(SET.Examples.caxis)
            hold on
            % Iterate over cluster and add boundaries
            for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
                % Get current cluster's pixels
                idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
                % Create a BW image of only the current cluster's
                % pixels being white
                BW = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
                for iPX = 1:length(idx)
                    BW(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 2)) = 1;
                end%iPX
                BW = BW';
                % Draw cluster's outlines
                [B, L] = bwboundaries(BW, 'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
                end%k
                clear k B L BW idx
            end%iCluster
            clear cm
        else
            caxis([min(min(DATA.ClusterResponse_all.mat)), max(max(DATA.ClusterResponse_all.mat))])
        end%if
        cb = colorbar;
        cb.Location = 'southoutside';
        title(SET.OdorNames{iOdor})
        set(gca, 'Fontsize', 3)
        
        
        % ---(ClusterOverview: Example traces)---------------------------------
        vec = [0:79]*0.2-4.8;
        if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
            % Iterate over all clusters
            for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
                % Iterate over stimuli
                % Preallocation
                TC = zeros(length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor})), length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).vec));
                % Iterate over repetitions
                for iRep = 1:length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor}))
                    % Get current data
                    temp = DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).SmoothResponse.(SET.OdorNames{iOdor}){iRep};
                    % Find current cluster's pixels
                    idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
                    % Iterate over frames
                    for iFrame = 1:length(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).vec)
                        % Preallocation
                        PX = zeros(1, length(idx));
                        % Get each pixel's activity at current frame
                        for iPx = 1:length(idx)
                            PX(iPx) =  temp(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPx), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPx), 2), iFrame);
                        end%iPx
                        % Get average activity of all pixels
                        TC(iRep, iFrame) = nanmean(PX);
                        
                    end%iFrame
                end%iRep
                % Plot
                subplot(7, 7, 14+iCluster);
                hold on
                plot(vec, nanmean(TC), 'Color', SET.Color(iOdor, :), 'LineWidth', 1)
                %                     plot(vec, nanmean(TC)-(std(TC)/sqrt(size(TC, 1))), 'Color', SET.Color(iOdor, :), 'LineWidth', 1)
                %                     plot(vec, nanmean(TC)+(std(TC)/sqrt(size(TC, 1))), 'Color', SET.Color(iOdor, :), 'LineWidth', 1)
                plot([vec(1), vec(end)], [0, 0], 'k')
                title(['Cluster: ', num2str(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster))])
                xlim([-1, 2])
                ylim(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).TC_ylim)
                set(gca, 'YTick', -0.03:0.005:0.03)
                set(gca, 'Fontsize', 3)
            end%iCluster
        end%if
        
        
        % ---(ClusterOverview: Cluster responses)--------------------------
        % Iterate over cluster
        for iCluster = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1)
            figure(hFig.ClusterOverview)
            subplot(7, 7, 21+iCluster)
            hold on
            temp_avg = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg);
            temp_sem = std(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg)./sqrt(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg,1));
            if temp_avg(iCluster)>=0
                plot([iOdor+0.25, iOdor+0.25], [temp_avg(iCluster), temp_avg(iCluster) + temp_sem(iCluster)], 'k')
                plot([iOdor, iOdor+0.5], [temp_avg(iCluster) + temp_sem(iCluster), temp_avg(iCluster) + temp_sem(iCluster)], 'k')
                rectangle('Position', [iOdor, 0, 0.5, temp_avg(iCluster)], 'FaceColor', SET.Color(iOdor, :), 'EdgeColor', 'none')
            elseif temp_avg(iCluster)<0
                plot([iOdor+0.25, iOdor+0.25], [temp_avg(iCluster), temp_avg(iCluster) - temp_sem(iCluster)], 'k')
                plot([iOdor, iOdor+0.5], [temp_avg(iCluster) - temp_sem(iCluster), temp_avg(iCluster) - temp_sem(iCluster)], 'k')
                rectangle('Position', [iOdor, temp_avg(iCluster), 0.5, -temp_avg(iCluster)], 'FaceColor', SET.Color(iOdor, :), 'EdgeColor', 'none')
            end
            title({['Cluster: ', num2str(iCluster)]; ['Size: ', num2str(length(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster))), ' pixels']})
            xlim([0.75, length(SET.OdorNames)+0.75])
            plot([0.75, length(SET.OdorNames)+0.75], [0 0], 'k')
            % Special treatment for examples
            if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
                ylim(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).ylim)
                set(gca, 'YTick', -0.03:0.005:0.03)
            else
                temp = [min(min(DATA.ClusterResponse_all.mat)), max(max(DATA.ClusterResponse_all.mat))];
                temp(1) = temp(1)-range(temp)*0.05;
                temp(2) = temp(2)+range(temp)*0.05;
                ylim(temp)
                clear temp
            end%if
            clear temp cm iPx ax
            set(gca, 'Fontsize', 3)
        end%iCluster
    end%iOdor
    
    
    % ---(ClusterCorrelation)----------------------------------------------
    figure(hFig.ClusterCorrelation)
    % --- overview
    ax(1,iFile) = subplot(length(SET.FileNames), 3, 1+(iFile-1)*3);
    imagesc(nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat,3))
    xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterMedoids, 1)+0.5])
    ylim([0.5, length(SET.OdorNames)+0.5])
    box off
    set(gca, 'XTick', [], 'YTick', [])
    colormap(ax(1,iFile), viridis(50))
    % Special treatment for examples
    if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
        caxis(SET.Examples.caxis)
    else
        caxis([min(min(DATA.ClusterResponse_all.mat)), max(max(DATA.ClusterResponse_all.mat))])
    end%if
    cb = colorbar;
    cb.Location = 'southoutside';
    % --- correlation
    ax(2,iFile) = subplot(length(SET.FileNames), 3, 2+(iFile-1)*3);
    imagesc(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.corrcoef)
    axis equal
    xlim([0.5, length(SET.OdorNames)+0.5])
    ylim([0.5, length(SET.OdorNames)+0.5])
    box off
    set(gca, 'XTick', [], 'YTick', [])
    colormap(ax(2,iFile), flipud(bone(50)))
    caxis([min(DATA.ClusterResponse_all.min_corrcoef) 1])
    cb = colorbar;
    cb.Location = 'southoutside';
    % --- dendrogram
    ax(3,iFile) = subplot(length(SET.FileNames), 3, 3+(iFile-1)*3);
    Y = pdist(nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.mat,3), 'euclidean');
    Z = linkage(Y, 'average');
    leafOrder = optimalleaforder(Z, Y);
    dendrogram(Z, 'Reorder', leafOrder)
    ylim([0 0.1])
    set(gca, 'XTickLabels', SET.OdorNames, 'view', [90 90])
    title(['c = ', num2str(cophenet(Z, Y))])
    
    
    % ---(save figures)----------------------------------------------------
    figure(hFig.PhysiologyHeatmaps)
    export_fig(['FIG\raw_imaging\PhysiologyHeatmaps_', SET.FileNames{iFile}(1:end-4)], '-pdf')
    figure(hFig.ClusterOverview)
    export_fig(['FIG\raw_imaging\ClusterOverview_', SET.FileNames{iFile}(1:end-4)], '-pdf')
    
end%iFile
clear iOdor iCluster ax*


% ---(ClusterCorrelation_all)----------------------------------------------
hFig.ClusterCorrelation_all = figure('Units', 'Normalized', 'Position', [0 0 1 1], 'Color', 'w');
% --- overview
ax(1) = subplot(1, 3, 1);
imagesc(DATA.ClusterResponse_all.mat)
xlim([0.5, size(DATA.ClusterResponse_all.mat, 2)+0.5])
ylim([0.5, size(DATA.ClusterResponse_all.mat, 1)+0.5])
box off
set(gca, 'XTick', [], 'YTick', [])
colormap(ax(1), viridis(50))
caxis([min(min(DATA.ClusterResponse_all.mat)), max(max(DATA.ClusterResponse_all.mat))])
cb = colorbar;
cb.Location = 'southoutside';
hold on
% Add separating lines
for iFile = 1:length(SET.FileNames)-1
    plot([SET.kmedoids.nCluster*iFile, SET.kmedoids.nCluster*iFile], [0.5, length(SET.OdorNames)+0.5], 'w', 'LineWidth', 2)
end
% --- correlation
ax(2) = subplot(1,3,2);
imagesc(DATA.ClusterResponse_all.corrcoef)
axis equal
xlim([0.5, length(SET.OdorNames)+0.5])
ylim([0.5, length(SET.OdorNames)+0.5])
box off
set(gca, 'XTick', [], 'YTick', [])
colormap(ax(2), flipud(bone(50)))
caxis([min(DATA.ClusterResponse_all.min_corrcoef) 1])
cb = colorbar;
cb.Location = 'southoutside';
% --- dendrogram
ax(3) = subplot(1,3,3);
Y = pdist(DATA.ClusterResponse_all.mat, 'euclidean');
Z = linkage(Y, 'average');
leafOrder = optimalleaforder(Z, Y);
dendrogram(Z, 'Reorder', leafOrder)
set(gca, 'XTickLabels', SET.OdorNames, 'view', [90 90])
title(['c = ', num2str(cophenet(Z, Y))])
clear ax


% ---(save figures)--------------------------------------------------------
figure(hFig.ClusterCorrelation)
export_fig(['FIG\raw_imaging\ClusterCorrelation'], '-pdf')
figure(hFig.ClusterCorrelation_all)
export_fig(['FIG\raw_imaging\ClusterCorrelation_all'], '-pdf')


% ---(ClusterOverview: Cluster responses)--------------------------
hFig.ExampleResponses = figure('Units', 'Normalized', 'Position', [0 0 1 1], 'Color', 'w');
for iFile = 1:2
    Val = [];
    subplot(2,1,iFile)
    hold on
    cnt=1;
    % Iterate over cluster
    for iCluster = SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster
        % Iterate over all stimuli
        for iOdor = 1:length(SET.OdorNames)
            
            temp_avg = nanmean(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(:,iCluster));
            temp_sem = std(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(:,iCluster))./sqrt(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterResponse.(SET.OdorNames{iOdor}).avg(:,iCluster),1));
            if temp_avg>=0
                plot([((cnt-1)*5)+iOdor+0.25, ((cnt-1)*5)+iOdor+0.25], [temp_avg, temp_avg + temp_sem], 'k')
                plot([((cnt-1)*5)+iOdor, ((cnt-1)*5)+iOdor+0.5], [temp_avg + temp_sem, temp_avg + temp_sem], 'k')
                rectangle('Position', [((cnt-1)*5)+iOdor, 0, 0.5, temp_avg], 'FaceColor', SET.Color(iOdor, :), 'EdgeColor', 'none')
            elseif temp_avg<0
                plot([((cnt-1)*5)+iOdor+0.25, ((cnt-1)*5)+iOdor+0.25], [temp_avg, temp_avg - temp_sem], 'k')
                plot([((cnt-1)*5)+iOdor, ((cnt-1)*5)+iOdor+0.5], [temp_avg - temp_sem, temp_avg - temp_sem], 'k')
                rectangle('Position', [((cnt-1)*5)+iOdor, temp_avg, 0.5, -temp_avg], 'FaceColor', SET.Color(iOdor, :), 'EdgeColor', 'none')
            end
            
            if temp_avg>0
                Val = [Val, temp_avg+temp_sem];
            else
                Val = [Val, temp_avg-temp_sem];
            end
            
            title({['Cluster: ', num2str(iCluster)]; ['Size: ', num2str(length(find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID==iCluster))), ' pixels']})
            
            clear temp cm iPx ax
        end%iCluster
        cnt = cnt+1;
    end%iOdor
    
    if min(Val) ==0
        %         ylim([0, max(Val)+range([min(Val), max(Val)])*0.01]);
    else
        %         ylim([min(Val)-range([min(Val), max(Val)])*0.01, max(Val)+range([min(Val), max(Val)])*0.01]);
    end
    
    xlim([0.75, ((cnt-2)*5)+iOdor+0.75])
    plot([0.75, ((cnt-2)*5)+iOdor+0.75], [0 0], 'k')
    
end

figure(hFig.ExampleResponses)
export_fig(['FIG\raw_imaging\ExampleResponses'], '-pdf')




%% Variability in cluster procedure

% Settings aquired post-clustering for examples
SET.Examples.file_181217_mp03_dff.Cluster = [1, 11, 3, 4, 24, 25];
SET.Examples.file_181219_mp04_dff.Cluster = [1, 3, 6, 12, 24, 25];
clc
numRep = 20;
clear ClusterID ClusterMedoids InterClusterDists ClusterVar

% Iterate over files
for iFile = 4:length(SET.FileNames)
    
    % Iterate over number of set replicates
    for iRep = 1:numRep
        [ClusterID{iRep,iFile}, ClusterMedoids{iRep,iFile}] = ...
            kmedoids(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.data_PCA, ...               data
            SET.kmedoids.nCluster, ...                                                                        clusters
            'OnlinePhase', SET.kmedoids.OnlinePhase, ...                                                      Turning on this feature tends to improve the quality of solutions
            'Distance', SET.kmedoids.Dist, ...                                                                'seuclidean' | 'cityblock' | 'minkowski' | 'chebychev' | 'mahalanobis' | 'cosine' | 'correlation' | 'spearman' | 'hamming' | 'jaccard'
            'Options', ...
            statset('MaxIter', SET.kmedoids.nIter, ...                                                        Maximum number of iterations allowed
            'UseParallel', SET.kmedoids.UseParallel), ...                                                     If true, compute in parallel
            'Replicates', SET.kmedoids.nRep, ...                                                              Number of times to repeat clustering using new initial cluster positions.
            'Start', SET.kmedoids.Start...                                                                   'plus' | 'sample' | 'cluster'
            );
        % Get the avergae between-cluster distance
        InterClusterDists(iRep,iFile) = nanmean(pdist(ClusterMedoids{iRep,iFile},'euclidean'));
        % Display current rep & file
        disp(['Done with rep. ', num2str(iRep), ' in file ', num2str(iFile)])
        
    end%iRep
    
    % Preallocation
    ClusterVar = ones(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF,1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF,2))*(-1/51);
    % Iterate over pixels
    for iPx = 1:size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates,1)
        % Get each replicate's medoid
        medoids = [];
        for iRep = 1:numRep
            medoids = [medoids;ClusterMedoids{iRep,iFile}(ClusterID{iRep,iFile}(iPx),:)];
        end%i
        % Get the average within-pixel distance between replicates
        ClusterVar(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(iPx,1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(iPx,2)) = nanmean(pdist(medoids,'euclidean'))/nanmean(InterClusterDists(:,iFile));
    end%iPx
    
    % Plot
    figure
    imagesc(ClusterVar')
    axis equal
    xlim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 2)+0.5])
    ylim([0.5, size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF', 1)+0.5])
    box off
    set(gca, 'XTick', [], 'YTick', [])
    cm = [0 0 0; viridis(length(unique(ClusterVar))*2)];
    colormap(cm)
    cb = colorbar;
    cb.Location = 'eastoutside';
    caxis([-1/51 1])
    cb.Ticks = [0 0.5 1];
    title(['file_', SET.FileNames{iFile}(1:end-4)], 'interpreter', 'none')
    % Add boundaries
    if strcmp(SET.FileNames{iFile}(1:end-4), '181217_mp03_dff') || strcmp(SET.FileNames{iFile}(1:end-4), '181219_mp04_dff')
        hold on
        % Iterate over cluster and add boundaries
        for iCluster = 1:length(SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster)
            % Get current cluster's pixels
            idx = find(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterID == SET.Examples.(['file_', SET.FileNames{iFile}(1:end-4)]).Cluster(iCluster));
            % Create a BW image of only the current cluster's
            % pixels being white
            BW = zeros(size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 1), size(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).BF, 2));
            for iPX = 1:length(idx)
                BW(DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 1), DATA.(['file_', SET.FileNames{iFile}(1:end-4)]).ClusterTabel.coordinates(idx(iPX), 2)) = 1;
            end%iPX
            BW = BW';
            % Draw cluster's outlines
            [B, L] = bwboundaries(BW, 'noholes');
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 1)
            end%k
            clear k B L BW idx
        end%iCluster
    end% if example
    
    % Save figure
    export_fig(['FIG\raw_imaging\ClusterVar_',SET.FileNames{iFile}(1:end-4)], '-pdf')
    
end%iFile









