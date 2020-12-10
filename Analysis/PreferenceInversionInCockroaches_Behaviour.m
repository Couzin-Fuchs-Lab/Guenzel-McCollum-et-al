%% Social modulation of individual preferences in cockroaches
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
% This script represents the pipeline used to analyse behavioural data.
% It first pools data across different trials, then analyses them in terms
% of, for example, how much time animals spent in any of the two shelters,
% last the script creates figures and calculates statistics.
%
% (MATLAB 2020a)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add path to save figures
addpath(genpath('altmany-export_fig'))


% Create folder to save plots
if ~isdir('FIG\raw'); mkdir('FIG\raw'); end


% Get ready
clear all
close all
clc


%%                     ----- SETTINGS -----
%**************************************************************************


% Information about videos
%--------------------------------------------------------------------------
% Set where the manual annotation of trials is located. This file should
% list the location and size of odor and control shelters, the center of
% the arena and its radius
SET.Path_Annotation = 'DATA\Annotations.xlsx';
% Import annotation
ANNOTATION = readtable(SET.Path_Annotation);
% Extract conditions names
SET.ConditionNames = unique(ANNOTATION.Condition);
% Extract odor names
SET.OdorNames = unique(ANNOTATION.Odor);
% Extract concentration names
SET.ConcentrationNames = unique(ANNOTATION.Concentration);

% Specify in which folder the tracked locations are saved
SET.Path_Data = 'Tracks\';

% Specify in which folder the tracked locations of only sheltered animals
% are saved (applies only for group trials)
SET.Path_Data_Shelter = 'BM3000\';

% Set appendix of files containing tracked locations
SET.FileName = '_compressed_tracked1.csv';

% Set appendix of files containing tracked locations of only sheltered
% animals (applies only for group trials)
SET.FileName_Shelter = '_shelter_tracked.csv';

% Give frame rate of recordings.
SET.FrameRate = 25;

% Set radius of shelters in (cm). This will be used to calculate the speed
% of the animal in cm/s instead of px/s
SET.ShelterRadius = 6; %[cm]



% Settings for data processing
% -------------------------------------------------------------------------
% Set how many minutes should be taken into account
SET.CutAfter = 30; %[min]
SET.CutAfter = SET.CutAfter*SET.FrameRate*60;

% Discard frames that are exceeding a certain distance to the arena's
% center. Note, the radius of the arena is 1.
SET.OutlierThreshold = 1.05;

% An animal might be right on the edge of the shelter. Set a tolerance to
% account for this by either increasing (>1) or decreasing (<1) the radii
% of the shelters
SET.DetectionTolerance = 7/6;

% Set minimum stay duration
SET.minStayDur = 0.5; %[sec]
SET.minStayDur = SET.minStayDur*SET.FrameRate;

% Set bins in which data will be divided for the heatmaps. Resulting number
% of bins will be SET.HeatMapGrid^2.
SET.HeatMapGrid = 100;

% Set number of bins for proportion of time
SET.PropTimeBins = 6;
SET.BinVec = linspace(0, SET.CutAfter, SET.PropTimeBins+1);
SET.BinVec(end) = SET.BinVec(end)+1;

% Set number of bins for insta entering freq
SET.EnterFreqTimeBins = 6;

% Cockroach motion is intermittent. Therefore, speed traces are filtered
% with a moving average.
% Set the width of the window used to smooth traces.
SET.SpeedWindow = 180; %[s]
SET.SpeedWindow = SET.SpeedWindow*SET.FrameRate;

% For some measures (e.g. heading), only propperly moving animals are
% taking into account
SET.minSpeed = 1; %[cm/s]
SET.maxSpeed = 50; %[cm/s]

% Set the bins used to calculate the relative distance to a shelter. Bins
% must have a range of [-1 1]
SET.HistBins = -1:0.1:1;

% Set number of resamples for bootstrapping.
SET.N_Boot = 10000;

% Test statistic for one-sample bootstrap hypothesis testing
SET.TestStat1 = @(x1, L_x1, PredetVal) abs(mean(x1) - PredetVal) ./ (std(x1) / L_x1);

% Test statistic for two-sample bootstrap hypothesis testing
SET.TestStat2 = @(x1, x2, L_x1, L_x2) abs(mean(x1) - mean(x2)) ./ sqrt(   std(x1)/L_x1 + std(x2)/L_x2   );

% Survival threshold
SET.SurvivalThreshold = 0.05;

% Settings for figure cosmetics
% -------------------------------------------------------------------------
% Set whether to plot in the first place
SET.DoPlot = 1;

% Set threshold to indicate border effects
SET.BorderOutlierThreshold = 0.85;

% Color code for different conditions
SET.Color = [...
    000   000   000;... GRP
    097   129   179;... Van3
    112   174   224;... Van2
    172   215   232;... Van2Fe3
    188   205   188;... Van2Fe2
    232   215   172;... Van3Fe2
    224   174   112;... Fe2
    179   129   097;... Fe3
    ]/255;

% Order of plots
SET.PlorOrder = {...
    'GRP VAN 0.01';...             GRP
    'IND VAN 0.001';...            Van3
    'IND VAN 0.01';...             Van2
    'IND FEC_VAN 0.001_0.01';...   Van2Fe3
    'IND FEC_VAN 0.01_0.01';...    Van2Fe2
    'IND FEC_VAN 0.01_0.001';...   Van3Fe2
    'IND FEC 0.01';...             Fe2
    'IND FEC 0.001';...            Fe3
    };



%%                         ----- GET DATA -----
%**************************************************************************


% Helper variables
cntTrial.Str = 'empty';
DATA.Duration = [];
maxSurvival = 0;
maxCount = [];
minCount = [];

% Iterate over all trials listed in the annotation
for iTrial = 1:size(ANNOTATION,1)
    
    % Get file
    curr.File = [ANNOTATION.Path{iTrial}, SET.Path_Data, ANNOTATION.Trial{iTrial}, SET.FileName];
    curr.Data = readtable(curr.File);
    
    % For tracking, videos have been scaled with a factor of 0.2. We account
    % for this with this variable.
    SET.ScalingFactor = ANNOTATION.ScaleFactor(iTrial);
    
    
    % Concentration
    curr.Concentration = ANNOTATION.Concentration{iTrial};
    curr.Concentration(curr.Concentration=='.') = '_';
    
    
    % Arena
    curr.Arena.center = eval(ANNOTATION.Arena_Center{iTrial});
    curr.Arena.radius = ANNOTATION.Arena_Radius(iTrial);
    
    
    % Shelter
    curr.Shelter.radius = ANNOTATION.Shelter_Radius(iTrial);
    curr.Shelter.Odor_Loc = eval(ANNOTATION.OdorShelter_Loc{iTrial});
    curr.Shelter.Control_Loc = eval(ANNOTATION.ControlShelter_Loc{iTrial});
    
    
    % Center data
    curr.Data.pos_x = curr.Data.pos_x*SET.ScalingFactor - curr.Arena.center(1);
    curr.Data.pos_y = curr.Data.pos_y*SET.ScalingFactor - curr.Arena.center(2);
    curr.Pos_all = [curr.Data.pos_x, curr.Data.pos_y];
    curr.Shelter.Odor_Loc(1) = curr.Shelter.Odor_Loc(1) - curr.Arena.center(1);
    curr.Shelter.Odor_Loc(2) = curr.Shelter.Odor_Loc(2) - curr.Arena.center(2);
    curr.Shelter.Control_Loc(1) = curr.Shelter.Control_Loc(1) - curr.Arena.center(1);
    curr.Shelter.Control_Loc(2) = curr.Shelter.Control_Loc(2) - curr.Arena.center(2);
    
    
    % Get number of individuals and divide data
    curr.IndNames = unique(curr.Data.id);
    curr.N = length(curr.IndNames);
    for iInd = 1:curr.N
        curr.Pos_X(:,iInd) = curr.Data.pos_x(find(strcmp(curr.Data.id, curr.IndNames{iInd})));
        curr.Pos_Y(:,iInd) = curr.Data.pos_y(find(strcmp(curr.Data.id, curr.IndNames{iInd})));
    end%iInd
    
    
    % Cut data to the same length (e.g. 25min)
    if SET.CutAfter > 0
        if length(curr.Pos_X) > SET.CutAfter
            curr.Pos_X(SET.CutAfter+1:end,:) = [];
            curr.Pos_Y(SET.CutAfter+1:end,:) = [];
        else
            shortby = SET.CutAfter-length(curr.Pos_X);
            curr.Pos_X = [curr.Pos_X; nan(shortby, curr.N)];
            curr.Pos_Y = [curr.Pos_Y; nan(shortby, curr.N)];
            clear shortby
        end
    end
    
    
    % Concatenate positional data of all individuals in one Mx2 matrix
    % with M = frames * number individuals
    curr.Pos_all = [reshape(curr.Pos_X, [size(curr.Pos_X,1)*size(curr.Pos_X,2), 1]), reshape(curr.Pos_Y, [size(curr.Pos_Y,1)*size(curr.Pos_Y,2), 1])];
    
    
    % Scale data to the radius of the arena
    curr.Pos_all = curr.Pos_all / ANNOTATION.Arena_Radius(iTrial);
    curr.Pos_X = curr.Pos_X / ANNOTATION.Arena_Radius(iTrial);
    curr.Pos_Y = curr.Pos_Y / ANNOTATION.Arena_Radius(iTrial);
    curr.Shelter.Odor_Loc = curr.Shelter.Odor_Loc / ANNOTATION.Arena_Radius(iTrial);
    curr.Shelter.Control_Loc = curr.Shelter.Control_Loc / ANNOTATION.Arena_Radius(iTrial);
    curr.Shelter.radius = curr.Shelter.radius / ANNOTATION.Arena_Radius(iTrial);
    
    
    % Exclude big outliers, i.e. where the tracking jumped.
    % --- For all data combined
    Dist2Center = sqrt(sum([curr.Pos_all(:, 1), curr.Pos_all(:, 2)]'.*[curr.Pos_all(:, 1), curr.Pos_all(:, 2)]'))';
    curr.Pos_all(Dist2Center > SET.OutlierThreshold, :) = NaN;
    % --- For each individual
    for iInd = 1:length(curr.IndNames)
        % Calculate distance to center
        Dist2Center = sqrt(sum([curr.Pos_X(:, iInd), curr.Pos_Y(:, iInd)]'.*[curr.Pos_X(:, iInd), curr.Pos_Y(:, iInd)]'))';
        % Discard outliers
        curr.Pos_X(Dist2Center > SET.OutlierThreshold, iInd) = NaN;
        curr.Pos_Y(Dist2Center > SET.OutlierThreshold, iInd) = NaN;
    end
    clear Dist2Center
    
    
    % Get the duration of the trial
    curr.Duration = size(curr.Pos_X, 1);
    
    
    % Rotate data
    % --- Get the current angle between the two shelters
    vec = curr.Shelter.Control_Loc - curr.Shelter.Odor_Loc;
    angle = atan2d(vec(2), vec(1));
    % --- Create rotation matrix
    M = [ cosd(-angle) -sind(-angle);
        sind(-angle)  cosd(-angle) ];
    % --- Rotate by multiplication with M
    curr.Shelter.Odor_Loc = (M*curr.Shelter.Odor_Loc(:,1:2)')';
    curr.Shelter.Control_Loc = (M*curr.Shelter.Control_Loc(:,1:2)')';
    curr.Pos_all = (M*curr.Pos_all(:,1:2)')';
    % Iterate over all individuals
    for iInd = 1:length(curr.IndNames)
        temp = (M*[curr.Pos_X(:,iInd), curr.Pos_Y(:,iInd)]')';
        curr.Pos_X(:,iInd) = temp(:,1);
        curr.Pos_Y(:,iInd) = temp(:,2);
        clear temp
    end%iInd
    clear vec angle M
    % Make sure that the odor location is right
    if curr.Shelter.Odor_Loc(1) < curr.Shelter.Control_Loc(1)
        angle = 180;
        % --- Create rotation matrix
        M = [ cosd(-angle) -sind(-angle);
            sind(-angle)  cosd(-angle) ];
        % --- Rotate by multiplication with M
        curr.Shelter.Odor_Loc = (M*curr.Shelter.Odor_Loc(:,1:2)')';
        curr.Shelter.Control_Loc = (M*curr.Shelter.Control_Loc(:,1:2)')';
        curr.Pos_all = (M*curr.Pos_all(:,1:2)')';
        % Iterate over all individuals
        for iInd = 1:length(curr.IndNames)
            temp = (M*[curr.Pos_X(:,iInd), curr.Pos_Y(:,iInd)]')';
            curr.Pos_X(:,iInd) = temp(:,1);
            curr.Pos_Y(:,iInd) = temp(:,2);
            clear temp
        end%iInd
        clear angle M
    end% if 180 rotation is needed
    
    
    % Pre-specification
    curr.Survival.Control_SurvivalOverTrial = [];
    curr.Survival.Odor_SurvivalOverTrial = [];
    curr.HeadingVector_all.x = [];
    curr.HeadingVector_all.y = [];
    curr.HeadingAngle_all = [];
    % Helper variables
    cnt.Survival.Odor_Exits = 1;
    cnt.Survival.Control_Exits = 1;
    cnt.Survival.Arena_Exits = 1;
    
    
    % For each individual, calculate (i) distance to shelters, (ii) whether
    % it was inside a shelter, (iii) the speed it was moving with, (iv)
    % heading angle
    for iInd = 1:length(curr.IndNames)
        
        
        % Caluclate angular location and save it as imaginary
        % number (angle+dist-i)
        theta = atan2d(curr.Pos_Y(:,iInd), curr.Pos_X(:,iInd));
        idx = find(theta<0);
        theta(idx) = 360-abs(theta(idx));
        helper = [curr.Pos_X(:,iInd), curr.Pos_Y(:,iInd)];
        dist = sqrt(sum(helper'.*helper'))';
        curr.PolarLoc(:,iInd) = complex(theta, dist);
        clear theta idx helper dist
        
        
        % Calculate the distance to the odor shelter
        helper = [curr.Pos_X(:,iInd), curr.Pos_Y(:,iInd)] - curr.Shelter.Odor_Loc;
        curr.Dist2Odor(:,iInd) = sqrt(sum(helper'.*helper'))';
        clear helper
        
        
        % Calculate the distance to the control shelter
        helper = [curr.Pos_X(:,iInd), curr.Pos_Y(:,iInd)] - curr.Shelter.Control_Loc;
        curr.Dist2Control(:,iInd) = sqrt(sum(helper'.*helper'))';
        clear helper
        
        
        % Check whether the animal was in one of the shelters
        curr.InShelter_Odor(:,iInd) = curr.Dist2Odor(:,iInd) < curr.Shelter.radius*SET.DetectionTolerance;
        curr.InShelter_Control(:,iInd) = curr.Dist2Control(:,iInd) < curr.Shelter.radius*SET.DetectionTolerance;
        
        
        % Exclude stay durations that are shorter than the required minimum
        % duration. Start with excluding the shortest possible stay of
        % [0,1,0], then exclude stays of increasing duration up to
        % SET.minStayDur.
        if SET.minStayDur > 0
            % --- Stay in odor shelter
            currStay = curr.InShelter_Odor(:,iInd);
            % Get periods of an animal being sheltered
            [L,num] = bwlabel(currStay);
            % Iterate over those periods
            for iStay = 1:num
                % Get location of period
                idx = find(L == iStay);
                % If the period, i.e. the stay, is too short, exclude it
                if length(idx) < SET.minStayDur
                    currStay(idx) = 0;
                end%if too short
            end%iStay
            curr.InShelter_Odor(:,iInd) = currStay;
            clear currStay idx iStay L num
            
            % --- Stay in control shelter
            currStay = curr.InShelter_Control(:,iInd);
            % Get periods of an animal being sheltered
            [L,num] = bwlabel(currStay);
            % Iterate over those periods
            for iStay = 1:num
                % Get location of period
                idx = find(L == iStay);
                % If the period, i.e. the stay, is too short, exclude it
                if length(idx) < SET.minStayDur
                    currStay(idx) = 0;
                end%if too short
            end%iStay
            curr.InShelter_Control(:,iInd) = currStay;
            clear currStay idx iStay L num
            
        end %if exclude short stays
        
        
        % Get shelter status for any shelter
        curr.InShelter_Any(:,iInd) = curr.InShelter_Odor(:,iInd) + curr.InShelter_Control(:,iInd);
        % Get unsheltered periods
        curr.InShelter_Arena(:,iInd) = curr.InShelter_Any(:,iInd) == 0;
        
        
        % Get how long animals 'survived' either in a shelter or outside of
        % it
        % --- Odor Exit
        idx_start=find(diff([0 curr.InShelter_Odor(:,iInd)' 0])==1);
        idx_end=find(diff([0 curr.InShelter_Odor(:,iInd)' 0])==-1);
        for iEnter = 1:length(idx_start)
            curr.Survival.Odor_Exits(cnt.Survival.Odor_Exits, 1:(idx_end(iEnter)-idx_start(iEnter))) = 1;
            cnt.Survival.Odor_Exits = cnt.Survival.Odor_Exits+1;
            
            % Also keep track whether survival time changes over time
            % of the trial. As we would like to save to variables at
            % ones, simply convert it to a complex number with the real
            % part representing the time of entry (min) and the
            % imaginary part the duration of stay (sec).
            curr.Survival.Odor_SurvivalOverTrial(end+1, 1) = complex(idx_start(iEnter)/(SET.FrameRate*60), (idx_end(iEnter)-idx_start(iEnter))/SET.FrameRate);
        end
        clear idx_start idx_end iEnter
        
        
        % --- Control Exit
        idx_start=find(diff([0 curr.InShelter_Control(:,iInd)' 0])==1);
        idx_end=find(diff([0 curr.InShelter_Control(:,iInd)' 0])==-1);
        for iEnter = 1:length(idx_start)
            curr.Survival.Control_Exits(cnt.Survival.Control_Exits, 1:(idx_end(iEnter)-idx_start(iEnter))) = 1;
            cnt.Survival.Control_Exits = cnt.Survival.Control_Exits+1;
            
            % Also keep track whether survival time changes over time
            % of the trial. As we would like to save to variables at
            % ones, simply convert it to a complex number with the real
            % part representing the time of entry (min) and the
            % imaginary part the duration of stay (sec).
            curr.Survival.Control_SurvivalOverTrial(end+1, 1) = complex(idx_start(iEnter)/(SET.FrameRate*60), (idx_end(iEnter)-idx_start(iEnter))/SET.FrameRate);
        end
        clear idx_start idx_end iEnter
        
        
        % --- Arena Exit
        idx_start=find(diff([0 curr.InShelter_Arena(:,iInd)' 0])==1);
        idx_end=find(diff([0 curr.InShelter_Arena(:,iInd)' 0])==-1);
        for iEnter = 1:length(idx_start)
            curr.Survival.Arena_Exits(cnt.Survival.Arena_Exits, 1:(idx_end(iEnter)-idx_start(iEnter))) = 1;
            cnt.Survival.Arena_Exits = cnt.Survival.Arena_Exits+1;
        end
        clear idx_start idx_end iEnter
        
        
        % Keep track of the size of the survival variable
        if isfield(curr.Survival, 'Odor_Exits') && size(curr.Survival.Odor_Exits,2) > maxSurvival
            maxSurvival = size(curr.Survival.Odor_Exits,2);
        end
        if isfield(curr.Survival, 'Control_Exits') && size(curr.Survival.Control_Exits,2) > maxSurvival
            maxSurvival = size(curr.Survival.Control_Exits,2);
        end
        if isfield(curr.Survival, 'Arena_Exits') && size(curr.Survival.Arena_Exits,2) > maxSurvival
            maxSurvival = size(curr.Survival.Arena_Exits,2);
        end
        
        
        % Combine the distance to shelters in one varible that is +1 if the
        % animal is in the center of the odor shelter, and -1 if it is in
        % the center of the contorl shelter
        % Ony keep distances for the closer shelter
        % --- Get indices
        idx_Odor =    find(curr.Dist2Odor(:,iInd) > curr.Dist2Control(:,iInd) | isnan(curr.Dist2Control(:,iInd)));
        helper.Odor = curr.Dist2Odor(:,iInd);
        idx_Control = find(curr.Dist2Odor(:,iInd) < curr.Dist2Control(:,iInd) | isnan(curr.Dist2Control(:,iInd)));
        helper.Control = curr.Dist2Control(:,iInd);
        % --- Remove
        helper.Odor(idx_Odor) = [];
        helper.Control(idx_Control) = [];
        % --- Scale to 1 (odor) and -1 (control)
        helper.Odor = abs((helper.Odor - max(helper.Odor)) / max(helper.Odor));
        helper.Control =  (helper.Control - max(helper.Control)) / max(helper.Control);
        % --- Combine
        curr.Dist2Both(:,iInd) = [helper.Odor; helper.Control; nan(sum(isnan(curr.Dist2Control(:,iInd))),1)];
        clear helper idx_Odor idx_Control
        
        
        % Calculate heading direction and speed
        helper_T2 = [curr.Pos_X(2:end,iInd), curr.Pos_Y(2:end,iInd)];
        helper_T1 = [curr.Pos_X(1:end-1,iInd), curr.Pos_Y(1:end-1,iInd)];
        helper = helper_T2 - helper_T1;
        
        
        % Get the speed, which is the length of the vetor between the
        % positions of two consecutive frames
        curr.Speed(:,iInd) = sqrt(sum(helper'.*helper'))';
        % If the animal was not moving
        curr.Speed(isnan(curr.Speed(:,iInd)),iInd) = 0;
        % Convert px/frame into cm/s
        curr.Speed(:,iInd) = ((curr.Speed(:,iInd)/curr.Shelter.radius)*SET.ShelterRadius)*SET.FrameRate;
        
        
        % Get the distance traveled
        curr.cumDistance(:,iInd) = cumsum(curr.Speed(:,iInd)/SET.FrameRate);
        
        
        % Get the heading as angle and vector
        curr.HeadingAngle(:,iInd) = atan2d(helper(:,2), helper(:,1));
        curr.HeadingAngle(curr.HeadingAngle(:,iInd)<0,iInd) = 180+abs(curr.HeadingAngle(curr.HeadingAngle(:,iInd)<0,iInd));
        curr.HeadingVector.x(:,iInd) = [helper(:,1); NaN];
        curr.HeadingVector.x(find(curr.Speed(:,iInd)<=SET.minSpeed | curr.Speed(:,iInd)>=SET.maxSpeed), iInd) = NaN;
        curr.HeadingVector.y(:,iInd) = [helper(:,2); NaN];
        curr.HeadingVector.y(find(curr.Speed(:,iInd)<=SET.minSpeed | curr.Speed(:,iInd)>=SET.maxSpeed), iInd) = NaN;
        clear helper_T0 helper_T1 helper
        
        
        % Pool direction for all individuals together in one variable, too
        curr.HeadingVector_all.x = [curr.HeadingVector_all.x; curr.HeadingVector.x(:,iInd)];
        curr.HeadingVector_all.y = [curr.HeadingVector_all.y; curr.HeadingVector.y(:,iInd)];
        curr.HeadingAngle_all = [curr.HeadingAngle_all; curr.HeadingAngle(:,iInd)];
        
    end%iInd
    clear cnt
    
    
    
    
    %----------------------------------------------------------------------
    % For group trials, sheltered individuals were re-tracked with constant
    % identities of entering and leaving animals. Note that those files
    % **only** contain the coordinates of sheltered animals.
    % Therefore, repeat parts of the script
    if strcmp(ANNOTATION.Condition{iTrial}, 'GRP')
        
        % Delete some of the previous calculations
        curr = rmfield(curr, 'InShelter_Odor');
        curr = rmfield(curr, 'InShelter_Control');
        curr = rmfield(curr, 'InShelter_Any');
        curr.Survival = rmfield(curr.Survival, 'Odor_Exits');
        curr.Survival = rmfield(curr.Survival, 'Control_Exits');
        curr.Survival.Control_SurvivalOverTrial = [];
        curr.Survival.Odor_SurvivalOverTrial = [];
        
        
        % Get file
        curr.File_Shelter = [ANNOTATION.Path{iTrial}, SET.Path_Data_Shelter, ANNOTATION.Trial{iTrial}, SET.FileName_Shelter];
        curr.Data_Shelter = readtable(curr.File_Shelter);
        
        
        % Get shelter locations
        curr.Shelter.Odor_Loc = eval(ANNOTATION.OdorShelter_Loc{iTrial});
        curr.Shelter.Control_Loc = eval(ANNOTATION.ControlShelter_Loc{iTrial});
        
        
        % Center data. No need to apply a scaling factor here
        curr.Data_Shelter.pos_x = curr.Data_Shelter.pos_x - curr.Arena.center(1);
        curr.Data_Shelter.pos_y = curr.Data_Shelter.pos_y - curr.Arena.center(2);
        curr.Shelter.Odor_Loc(1) = curr.Shelter.Odor_Loc(1) - curr.Arena.center(1);
        curr.Shelter.Odor_Loc(2) = curr.Shelter.Odor_Loc(2) - curr.Arena.center(2);
        curr.Shelter.Control_Loc(1) = curr.Shelter.Control_Loc(1) - curr.Arena.center(1);
        curr.Shelter.Control_Loc(2) = curr.Shelter.Control_Loc(2) - curr.Arena.center(2);
        
        % Get number of individuals and divide data
        curr.IndNames_Shelter = unique(curr.Data_Shelter.id);
        for iInd = 1:length(curr.IndNames_Shelter)
            curr.Pos_X_Shelter(:,iInd) = curr.Data_Shelter.pos_x(find(strcmp(curr.Data_Shelter.id, curr.IndNames_Shelter{iInd})));
            curr.Pos_Y_Shelter(:,iInd) = curr.Data_Shelter.pos_y(find(strcmp(curr.Data_Shelter.id, curr.IndNames_Shelter{iInd})));
        end%iInd
        
        
        % Cut data to the same length (e.g. 25min)
        if SET.CutAfter > 0
            if length(curr.Pos_X_Shelter) > SET.CutAfter
                curr.Pos_X_Shelter(SET.CutAfter+1:end,:) = [];
                curr.Pos_Y_Shelter(SET.CutAfter+1:end,:) = [];
            else
                shortby = SET.CutAfter-length(curr.Pos_X);
                curr.Pos_X_Shelter = [curr.Pos_X_Shelter; nan(shortby, curr.N)];
                curr.Pos_Y_Shelter = [curr.Pos_Y_Shelter; nan(shortby, curr.N)];
                clear shortby
            end
        end
        
        
        % Scale data to the radius of the arena
        curr.Pos_X_Shelter = curr.Pos_X_Shelter / ANNOTATION.Arena_Radius(iTrial);
        curr.Pos_Y_Shelter = curr.Pos_Y_Shelter / ANNOTATION.Arena_Radius(iTrial);
        curr.Shelter.Odor_Loc = curr.Shelter.Odor_Loc / ANNOTATION.Arena_Radius(iTrial);
        curr.Shelter.Control_Loc = curr.Shelter.Control_Loc / ANNOTATION.Arena_Radius(iTrial);
        
        
        % Rotate data
        % --- Get the current angle between the two shelters
        vec = curr.Shelter.Control_Loc - curr.Shelter.Odor_Loc;
        angle = atan2d(vec(2), vec(1));
        % --- Create rotation matrix
        M = [ cosd(-angle) -sind(-angle);
            sind(-angle)  cosd(-angle) ];
        % Iterate over all individuals
        for iInd = 1:length(curr.IndNames_Shelter)
            temp = (M*[curr.Pos_X_Shelter(:,iInd), curr.Pos_Y_Shelter(:,iInd)]')';
            curr.Pos_X_Shelter(:,iInd) = temp(:,1);
            curr.Pos_Y_Shelter(:,iInd) = temp(:,2);
            clear temp
        end%iInd
        % Rotate shelter locations
        curr.Shelter.Odor_Loc =    (M*curr.Shelter.Odor_Loc')';
        curr.Shelter.Control_Loc = (M*curr.Shelter.Control_Loc')';
        clear vec angle M
        
        % Make sure that the odor location is right
        if curr.Shelter.Odor_Loc(1) < curr.Shelter.Control_Loc(1)
            angle = 180;
            % --- Create rotation matrix
            M = [ cosd(-angle) -sind(-angle);
                sind(-angle)  cosd(-angle) ];
            % --- Rotate by multiplication with M
            curr.Shelter.Odor_Loc = (M*curr.Shelter.Odor_Loc(:,1:2)')';
            curr.Shelter.Control_Loc = (M*curr.Shelter.Control_Loc(:,1:2)')';
            % Iterate over all individuals
            for iInd = 1:length(curr.IndNames_Shelter)
                temp = (M*[curr.Pos_X_Shelter(:,iInd), curr.Pos_Y_Shelter(:,iInd)]')';
                curr.Pos_X_Shelter(:,iInd) = temp(:,1);
                curr.Pos_Y_Shelter(:,iInd) = temp(:,2);
                clear temp
            end%iInd
            clear angle M
        end% if 180 rotation is needed
        
        
        % Helper variables
        cnt.Survival.Odor_Exits = 1;
        cnt.Survival.Control_Exits = 1;
        cnt.Survival.Arena_Exits = 1;
        
        % For each individual, calculate whether it was inside a shelter
        for iInd = 1:length(curr.IndNames_Shelter)
            
            
            % Calculate the distance to the odor shelter
            helper = [curr.Pos_X_Shelter(:,iInd), curr.Pos_Y_Shelter(:,iInd)] - curr.Shelter.Odor_Loc;
            curr.Dist2Odor_Shelter(:,iInd) = sqrt(sum(helper'.*helper'))';
            clear helper
            
            
            % Calculate the distance to the control shelter
            helper = [curr.Pos_X_Shelter(:,iInd), curr.Pos_Y_Shelter(:,iInd)] - curr.Shelter.Control_Loc;
            curr.Dist2Control_Shelter(:,iInd) = sqrt(sum(helper'.*helper'))';
            clear helper
            
            
            % Check whether the animal was in one of the shelters
            curr.InShelter_Odor(:,iInd) = curr.Dist2Odor_Shelter(:,iInd) < curr.Shelter.radius*SET.DetectionTolerance;
            curr.InShelter_Control(:,iInd) = curr.Dist2Control_Shelter(:,iInd) < curr.Shelter.radius*SET.DetectionTolerance;
            
            
            % Exclude stay durations that are shorter than the required minimum
            % duration. Start with excluding the shortest possible stay of
            % [0,1,0], then exclude stays of increasing duration up to
            % SET.minStayDur.
            if SET.minStayDur > 0
                % --- Stay in odor shelter
                currStay = curr.InShelter_Odor(:,iInd);
                % Get periods of an animal being sheltered
                [L,num] = bwlabel(currStay);
                % Iterate over those periods
                for iStay = 1:num
                    % Get location of period
                    idx = find(L == iStay);
                    % If the period, i.e. the stay, is too short, exclude it
                    if length(idx) < SET.minStayDur
                        currStay(idx) = 0;
                    end%if too short
                end%iStay
                curr.InShelter_Odor(:,iInd) = currStay;
                clear currStay idx iStay L num
                
                % --- Stay in control shelter
                currStay = curr.InShelter_Control(:,iInd);
                % Get periods of an animal being sheltered
                [L,num] = bwlabel(currStay);
                % Iterate over those periods
                for iStay = 1:num
                    % Get location of period
                    idx = find(L == iStay);
                    % If the period, i.e. the stay, is too short, exclude it
                    if length(idx) < SET.minStayDur
                        currStay(idx) = 0;
                    end%if too short
                end%iStay
                curr.InShelter_Control(:,iInd) = currStay;
                clear currStay idx iStay L num
                
            end
            
            
            % Get shelter status for any shelter
            curr.InShelter_Any(:,iInd) = curr.InShelter_Odor(:,iInd) + curr.InShelter_Control(:,iInd);
            
            
            % Get how long animals 'survived' either in a shelter or outside of
            % it
            % --- Odor Exit
            idx_start=find(diff([0 curr.InShelter_Odor(:,iInd)' 0])==1);
            idx_end=find(diff([0 curr.InShelter_Odor(:,iInd)' 0])==-1);
            for iEnter = 1:length(idx_start)
                curr.Survival.Odor_Exits(cnt.Survival.Odor_Exits, 1:(idx_end(iEnter)-idx_start(iEnter))) = 1;
                cnt.Survival.Odor_Exits = cnt.Survival.Odor_Exits+1;
                
                % Also keep track whether survival time changes over time
                % of the trial. As we would like to save to variables at
                % ones, simply convert it to a complex number with the real
                % part representing the time of entry (min) and the
                % imaginary part the duration of stay (sec).
                curr.Survival.Odor_SurvivalOverTrial(end+1, 1) = complex(idx_end(iEnter)/(SET.FrameRate*60), (idx_end(iEnter)-idx_start(iEnter))/(SET.FrameRate*60));
            end
            clear idx_start idx_end iEnter
            
            
            % --- Control Exit
            idx_start=find(diff([0 curr.InShelter_Control(:,iInd)' 0])==1);
            idx_end=find(diff([0 curr.InShelter_Control(:,iInd)' 0])==-1);
            for iEnter = 1:length(idx_start)
                curr.Survival.Control_Exits(cnt.Survival.Control_Exits, 1:(idx_end(iEnter)-idx_start(iEnter))) = 1;
                cnt.Survival.Control_Exits = cnt.Survival.Control_Exits+1;
                
                % Also keep track whether survival time changes over time
                % of the trial. As we would like to save to variables at
                % ones, simply convert it to a complex number with the real
                % part representing the time of entry (min) and the
                % imaginary part the duration of stay (sec).
                curr.Survival.Control_SurvivalOverTrial(end+1, 1) = complex(idx_end(iEnter)/(SET.FrameRate*60), (idx_end(iEnter)-idx_start(iEnter))/(SET.FrameRate*60));
            end
            clear idx_start idx_end iEnter
            
            
            % Keep track of the size of the survival variable
            if isfield(curr.Survival, 'Odor_Exits') && size(curr.Survival.Odor_Exits,2) > maxSurvival
                maxSurvival = size(curr.Survival.Odor_Exits,2);
            end
            if isfield(curr.Survival, 'Control_Exits') && size(curr.Survival.Control_Exits,2) > maxSurvival
                maxSurvival = size(curr.Survival.Control_Exits,2);
            end
            
        end%iInd
        clear cnt
    end
    %----------------------------------------------------------------------
    
    
    
    
    % Counter
    if ~strcmp(cntTrial.Str, [ANNOTATION.Condition{iTrial},ANNOTATION.Odor{iTrial},['C', curr.Concentration]])
        % Pre-specification
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Pos_all = [];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Direction_all = [];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.x = [];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.y = [];
        % Counter
        cntTrial.Str = [ANNOTATION.Condition{iTrial},ANNOTATION.Odor{iTrial},['C', curr.Concentration]];
        cntTrial.Val = 1;
    end
    
    
    % Some trials were divided into two files. First trial name ends with
    % *T1 and the second with *T2. Make sure to concatenate those.
    curr.NewTrial = 1;
    if iTrial > 1 % very first entry in list can't be the second trial
        % Check whether current trial ID matched the previous
        if strcmp(ANNOTATION.Trial{iTrial-1}(1:end-2), ANNOTATION.Trial{iTrial}(1:end-2))
            % Check whether the current entry is a succeeding part of the
            % previous entry
            if strcmp(ANNOTATION.Trial{iTrial-1}(end-1:end), 'T1') && strcmp(ANNOTATION.Trial{iTrial}(end-1:end), 'T2')
                % Indicate that this is not a new trial
                curr.NewTrial = 0;
            end%if no new trial
        end%if matching IDs
    end%if not first entry in list
    
    
    % Sort data in one variable "DATA"
    if curr.NewTrial
        % --- Trial-wise
        DATA.Duration = [DATA.Duration; size(curr.Pos_X,1)];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).N =            	        curr.N;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Duration =            	curr.Duration;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_all =                 curr.Pos_all;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).PolarLoc =                curr.PolarLoc;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_X =                   curr.Pos_X;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_Y =                   curr.Pos_Y;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).cumDistance =             curr.cumDistance;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Speed =                   curr.Speed;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingAngle =            curr.HeadingAngle;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.x =         curr.HeadingVector.x;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.y =         curr.HeadingVector.y;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Shelter_Radius =          curr.Shelter.radius;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Odor_Loc =                curr.Shelter.Odor_Loc;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Control_Loc =             curr.Shelter.Control_Loc;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Odor =               curr.Dist2Odor;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Control =            curr.Dist2Control;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Both =               curr.Dist2Both;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Odor =          curr.InShelter_Odor;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Control =       curr.InShelter_Control;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Any =           curr.InShelter_Any;
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Survival =                curr.Survival;
        
        % --- Condition-wise
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).ID{cntTrial.Val,1} = ANNOTATION.Trial{iTrial};
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Pos_all = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Pos_all; ...
            curr.Pos_all];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Direction_all = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Direction_all;...
            curr.HeadingAngle_all];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.x = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.x;...
            curr.HeadingVector_all.x];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.y = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.y;...
            curr.HeadingVector_all.y];
        
        % Update counter
        cntTrial.Val = cntTrial.Val+1;
        
    else
        % Count back
        cntTrial.Val = cntTrial.Val-1;
        
        % Pool survival curves
        temp_Survival = DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Survival;
        % --- Odor_Exits
        for iSurv = 1:size(curr.Survival.Odor_Exits)
            temp_Survival.Odor_Exits(end+1, 1:length(curr.Survival.Odor_Exits(iSurv,:))) = curr.Survival.Odor_Exits(iSurv,:);
        end%isurv
        % --- Odor_Exits
        for iSurv = 1:size(curr.Survival.Control_Exits)
            temp_Survival.Control_Exits(end+1, 1:length(curr.Survival.Control_Exits(iSurv,:))) = curr.Survival.Control_Exits(iSurv,:);
        end%isurv
        % --- Arena_Exits
        for iSurv = 1:size(curr.Survival.Arena_Exits)
            temp_Survival.Arena_Exits(end+1, 1:length(curr.Survival.Arena_Exits(iSurv,:))) = curr.Survival.Arena_Exits(iSurv,:);
        end%isurv
        
        % --- Trial-wise
        DATA.Duration(end) = DATA.Duration(end) + size(curr.Pos_X,1);
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Duration =                DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Duration + size(curr.Pos_X,1);
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_all =                 [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_all;                     curr.Pos_all];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).PolarLoc =                [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).PolarLoc;                    curr.PolarLoc];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_X =                   [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_X;                       curr.Pos_X];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_Y =                   [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Pos_Y;                       curr.Pos_Y];
        for iInd = 1:size(DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).cumDistance,2)
            curr.cumDistance(:,iInd) = curr.cumDistance(:,iInd) + DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).cumDistance(end,iInd);
        end
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).cumDistance =             [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).cumDistance;                 curr.cumDistance];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Speed =                   [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Speed;                       curr.Speed];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingAngle =            [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingAngle;                curr.HeadingAngle];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.x =         [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.x;             curr.HeadingVector.x];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.y =         [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).HeadingVector.x;             curr.HeadingVector.y];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Shelter_Radius =          nanmean([DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Shelter_Radius;      curr.Shelter.radius]);
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Odor_Loc =                nanmean([DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Odor_Loc;            curr.Shelter.Odor_Loc]);
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Control_Loc =             nanmean([DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Control_Loc;         curr.Shelter.Control_Loc]);
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Odor =               [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Odor;                   curr.Dist2Odor];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Control =            [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Control;                curr.Dist2Control];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Both =               [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Dist2Both;                   curr.Dist2Both];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Odor =          [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Odor;              curr.InShelter_Odor];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Control =       [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Control;           curr.InShelter_Control];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Any =           [DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).InShelter_Any;               curr.InShelter_Any];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).(['Trial_', num2str(cntTrial.Val)]).Survival =                temp_Survival; clear temp_Survival
        
        % --- Condition-wise
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).ID{cntTrial.Val,1} = ANNOTATION.Trial{iTrial};
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Pos_all = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Pos_all;...
            curr.Pos_all];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Direction_all = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).Direction_all;...
            curr.HeadingAngle_all];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.x = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.x;...
            curr.HeadingVector_all.x];
        DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.y = [...
            DATA.(ANNOTATION.Condition{iTrial}).(ANNOTATION.Odor{iTrial}).(['C', curr.Concentration]).HeadingVector_all.y;...
            curr.HeadingVector_all.y];
        
        % Update counter
        cntTrial.Val = cntTrial.Val+1;
        
    end%if new trial
    
    clear curr
    
end%iTrial

% Before continuing with pooling the data, save the current workspace if
% something breaks later
save('tempDATA.mat', 'DATA', 'ANNOTATION', 'SET', 'maxSurvival', '-v7.3')

% Delete what has been forgotton to delete on the way
clear tExits tEntries nExits nEntries iTrial iSurv iInd iEntry helper_T2 Exits Entries d_maxSurvival d_maxDur cntTrial cnt pattval iStart iPattern

%%                       ----- POOL DATA -----
%**************************************************************************


% Get the duration of the longest trial to have equal x axis
maxDur = max(DATA.Duration);

counterTable = 1;

% Iterate over all conditions
for cntCond = 1:length(SET.PlorOrder)
    
    
    % Get the right combination of condition, odor, and concentration
    currPlot = SET.PlorOrder(cntCond,:);
    idx = cell2mat(strfind(currPlot, ' '));
    iCondition = find(strcmp(SET.ConditionNames, currPlot{1}(1:idx(1)-1)));
    iOdor = find(strcmp(SET.OdorNames, currPlot{1}(idx(1)+1:idx(2)-1)));
    iConcentration = find(strcmp(SET.ConcentrationNames, currPlot{1}(idx(2)+1:end)));
    clear idx currPlot
    
    
    % Create a string to identify the current concentration
    Concentration = ['C', SET.ConcentrationNames{iConcentration}];
    Concentration(Concentration=='.') = '_';
    
    
    % Continue, if data of the current condition, odor, and concentration have been pooled previously
    if ~isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}), Concentration)
        error('Combination of condition, odor, and concentration does not exist')
    end
    
    
    % Avg over trials. Therefore, prepecify variables
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Dist2Both = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitIndex = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed_movOnly = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Delta = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Any = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Delta = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Any = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).OdorPos = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ControlPos = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterRadius = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InstaVisitFreq_Any = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InstaVisitFreq_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InstaVisitFreq_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InstaVisitFreq_Delta = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Any = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Odor = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Control = [];
    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Delta = [];
    
    
    % Introduce several counter variables
    clear cnt
    cnt.Odor_Entry = 1;
    cnt.Control_Entry = 1;
    cnt.Odor_Exit = 1;
    cnt.Control_Exit = 1;
    
    
    % Iterate over all trials of the current
    % combination of condition+odor+concentration
    for iTrial = 1:length(fields(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration)))-4
        
        d_maxDur = maxDur - size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Dist2Both, 1);
        
        if strcmp(SET.ConditionNames{iCondition}, 'GRP')
            
            % --- I-Entry-I (Odor)
            Entries.Odor = diff(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,2));
            Entries.Odor(Entries.Odor<0)=0;
            tEntries.Odor = find(Entries.Odor>0);
            nEntries.Odor = length(tEntries.Odor);
            for iEntry = 2:nEntries.Odor
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Odor(1,cnt.Odor_Entry) = tEntries.Odor(iEntry) - tEntries.Odor(iEntry-1);
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Odor(2,cnt.Odor_Entry) = sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor(tEntries.Odor(iEntry)-1,:));
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Odor(3,cnt.Odor_Entry) = iTrial;
                cnt.Odor_Entry = cnt.Odor_Entry+1;
            end%iEntry
            
            % --- I-Entry-I (Control)
            Entries.Control = diff(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,2));
            Entries.Control(Entries.Control<0)=0;
            tEntries.Control = find(Entries.Control>0);
            nEntries.Control = length(tEntries.Control);
            for iEntry = 2:nEntries.Control
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Control(1,cnt.Control_Entry) = tEntries.Control(iEntry) - tEntries.Control(iEntry-1);
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Control(2,cnt.Control_Entry) = sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control(tEntries.Control(iEntry)-1,:));
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterEntryInterval.Control(3,cnt.Control_Entry) = iTrial;
                cnt.Control_Entry = cnt.Control_Entry+1;
            end%iEntry
            
            % --- I-Exit-I (Odor)
            Exits.Odor = diff(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,2));
            Exits.Odor(Exits.Odor>0)=0;
            tExits.Odor = find(Exits.Odor<0);
            nExits.Odor = length(tExits.Odor);
            for iEntry = 2:nExits.Odor
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Odor(1,cnt.Odor_Exit) = tExits.Odor(iEntry) - tExits.Odor(iEntry-1);
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Odor(2,cnt.Odor_Exit) = sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor(tExits.Odor(iEntry)-1,:));
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Odor(3,cnt.Odor_Exit) = iTrial;
                cnt.Odor_Exit = cnt.Odor_Exit+1;
            end%iEntry
            
            % --- I-Exit-I (Control)
            Exits.Control = diff(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,2));
            Exits.Control(Exits.Control>0)=0;
            tExits.Control = find(Exits.Control<0);
            nExits.Control = length(tExits.Control);
            for iEntry = 2:nExits.Control
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Control(1,cnt.Control_Exit) = tExits.Control(iEntry) - tExits.Control(iEntry-1);
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Control(2,cnt.Control_Exit) = sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control(tExits.Control(iEntry)-1,:));
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).InterExitInterval.Control(3,cnt.Control_Exit) = iTrial;
                cnt.Control_Exit = cnt.Control_Exit+1;
            end%iEntry
            
        end%if grp trial
        
        
        % --- Dist2Both
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Dist2Both = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Dist2Both;...
            hist(reshape(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Dist2Both, [numel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Dist2Both), 1]), SET.HistBins) / numel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Dist2Both)];
        
        
        % --- StayDur_Odor
        [L,num] = bwlabel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor);
        if num>0
            dur = nan(1,num);
            for iStay = 1:num
                dur(iStay) = sum(sum(L==iStay))/SET.FrameRate;
            end
            temp = mean(dur);
        else
            temp = 0;
        end
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor;...
            temp];
        clear L num dur iStay temp
        
        
        % --- StayDur_Control
        [L,num] = bwlabel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control);
        if num>0
            dur = nan(1,num);
            for iStay = 1:num
                dur(iStay) = sum(sum(L==iStay))/SET.FrameRate;
            end
            temp = mean(dur);
        else
            temp = 0;
        end
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control;...
            temp];
        clear L num dur iStay temp
        
        
        % --- EntryShelter_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor;...
            sum(sum(diff(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,1,1)==1))];
        
        
        % --- EntryShelter_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control...
            sum(sum(diff(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,1,1)==1))];
        
        
        % --- ExitShelter_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor;...
            sum(sum(diff(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,1,1)==-1))];
        
        
        % --- ExitShelter_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control;...
            sum(sum(diff(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,1,1)==-1))];
        
        
        % --- EntryIndex
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex;...
            (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(end)-PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(end)) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(end))];
        
        
        % --- ExitIndex
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitIndex = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitIndex;...
            (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(end)-PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(end)) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(end))];
        
        
        % --- PropEntryShelter_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Odor;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(end) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(end))];
        
        
        % --- PropEntryShelter_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Control;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(end) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(end))];
        
        
        % --- PropExitShelter_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Odor;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(end) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(end))];
        
        
        % --- PropExitShelter_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Control;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(end) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(end))];
        
        
        % --- Speed
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed;...
            mean(movmean([DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Speed; nan(d_maxDur,1)]', SET.SpeedWindow, 2, 'omitnan'),1)];
        
        
        % --- Speed_movOnly
        temp = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Speed;
        temp(temp==0)=NaN;
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed_movOnly = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed_movOnly;...
            mean(movmean([temp; nan(d_maxDur,1)]', SET.SpeedWindow, 2, 'omitnan'),1)];

        
        % --- PropTimeSheltered_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor;...
            sum(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor)) / numel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor)];
        
        
        % --- PropTimeSheltered_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control;...
            sum(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control)) / numel(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control)];
        
        
        % --- ShelterTimeIndex
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex;...
            (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor(end) - PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control(end)) / (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor(end) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control(end))];

        
        % --- cumPropTimeSheltered_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor;...
            cumsum(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,2)  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,2))'  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor,1)];
        
        
        % --- cumPropTimeSheltered_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control;...
            cumsum(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,2)  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,2))'  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control,1)];
        
        
        % --- cumPropTimeSheltered_Delta
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Delta = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Delta;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor(end,:) - PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control(end,:)];
        
        
        % --- cumPropTimeSheltered_Any
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Any = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Any;...
            cumsum(sum(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Any,2)  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Any,2))'  /  size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Any,1)'];
        
        
        % --- normcumPropTimeSheltered_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Odor;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor(end, :) / PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Any(end, end)];
        
        
        % --- normcumPropTimeSheltered_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Control;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control(end, :) / PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Any(end, end)];
        
        
        % --- normcumPropTimeSheltered_Delta
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Delta = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Delta;...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Odor(end,:) - PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Control(end,:)];
        
        
        % --- binnedPropTime_Any
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        AnyData = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Any;
        % Check for more than one animals
        if size(AnyData,2)>1
            AnyData = sum(AnyData,2)/size(AnyData,2);
        end
        AnyData = cumsum(AnyData)/size(AnyData,1);
        % Bin data into n bins
        [~,idx] = histc(vec, SET.BinVec);
        AnyData_bin = accumarray(idx(:), AnyData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Any = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Any; ...
            AnyData_bin'];
        
        
        % --- binnedPropTime_Odor
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        OdorData = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Odor;
        if size(OdorData,2)>1
            OdorData = sum(OdorData,2)/size(OdorData,2);
        end
        OdorData = cumsum(OdorData)/size(OdorData,1);
        % Bin data into X bins
        [~,idx] = histc(vec, SET.BinVec);
        OdorData_bin = accumarray(idx(:), OdorData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Odor; ...
            OdorData_bin'];
        
        
        % --- binnedPropTime_Control
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        ControlData = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).InShelter_Control;
        if size(ControlData,2)>1
            ControlData = sum(ControlData,2)/size(ControlData,2);
        end
        ControlData = cumsum(ControlData)/size(ControlData,1);
        % Bin data into X bins
        [~,idx] = histc(vec, SET.BinVec);
        ControlData_bin = accumarray(idx(:), ControlData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Control; ...
            ControlData_bin'];
        
        
        % --- binnedPropTime_Delta
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta; ...
            OdorData_bin'-ControlData_bin'];
        clear vec OdorData* ControlData* idx
             
        
        % --- Survival_Odor
        if isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival, 'Odor_Exits')
            d_maxSurvival = maxSurvival - size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits, 2);
            if size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits, 1) > 1
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor;...
                    nanmean([DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits, 1), d_maxSurvival)])];
            else
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor;...
                    [DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_Exits, 1), d_maxSurvival)]];
            end%if more than one exit
        else
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor;...
                nan(1,SET.CutAfter)];
        end%if Odor_Exits exists
        
        
        % --- Survival_Control
        if isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival, 'Control_Exits')
            d_maxSurvival = maxSurvival - size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits, 2);
            if size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits, 1) > 1
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control;...
                    nanmean([DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits, 1), d_maxSurvival)])];
            else
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control;...
                    [DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_Exits, 1), d_maxSurvival)]];
            end%if more than one exit
        else
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control;...
                nan(1,SET.CutAfter)];
        end%if Control_Exits exists
        
        
        % --- Survival_Arena
        if isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival, 'Arena_Exits')
            d_maxSurvival = maxSurvival - size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits, 2);
            if size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits, 1) > 1
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena;...
                    nanmean([DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits, 1), d_maxSurvival)])];
            else
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena = [...
                    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena;...
                    [DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits,...
                    zeros(size(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Arena_Exits, 1), d_maxSurvival)]];
            end%if more than one exit
        else
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Arena;...
                nan(1,SET.CutAfter)];
        end%if Arena_Exits exists
        
        
        % --- SurvivalThreshold_Odor
        if isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival, 'Odor_Exits')
            currSurv = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(iTrial,:);
            [~, idx] = min(abs(currSurv - SET.SurvivalThreshold));
        else
            idx = NaN;
        end
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Odor;...
            idx /(SET.FrameRate*60)];
        clear currSurv idx
        
        
        % --- SurvivalThreshold_Control
        if isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival, 'Control_Exits')
            currSurv = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(iTrial,:);
            [~, idx] = min(abs(currSurv - SET.SurvivalThreshold));
        else
            idx = NaN;
        end
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Control;...
            idx /(SET.FrameRate*60)]; 
        clear currSurv idx
        
        
        % --- SurvivalTimeOverTrial_Odor
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Odor;...
            DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Odor_SurvivalOverTrial];
        
        
        % --- SurvivalTimeOverTrial_Control
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalTimeOverTrial_Control;...
            DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Survival.Control_SurvivalOverTrial];
        
        
        % --- SurvivalIndex
        if isnan(sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(end,:)))
            
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex;...
                (0 - sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(end,:))) / (0 + sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(end,:)))];
            
        elseif isnan(sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(end,:)))
            
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex;...
                (sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(end,:)) - 0) / (sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(end,:)) + 0)];
            
        else
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex = [...
                PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex;...
                (sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(end,:)) - sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(end,:))) / (sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor(end,:)) + sum(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control(end,:)))];
        end
        
        
        % --- Location of odor shelter
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).OdorPos = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).OdorPos;...
            DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Odor_Loc];
        
        
        % --- Location of control shelter
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ControlPos = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ControlPos;...
            DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Control_Loc];
        
        
        % --- Radius of shelters
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterRadius = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterRadius;...
            DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Shelter_Radius];
        
        
        % --- All positions
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
        
        
        % --- All headings
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).HeadingVector_all = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).HeadingVector_all;
        
        
        % --- Prefered location in arena
        % Get current position
        currPos = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_',num2str(iTrial)]).Pos_all;
        % Get location as angle
        currTheta = atan2d(currPos(:,2), currPos(:,1));
        currTheta(currTheta<0) = 360-abs(currTheta(currTheta<0));
        % Get position to centre
        Dist2Centre = sqrt(sum(currPos'.*currPos'))';
        % Define several logicals to later sort
        OK.general = ~isnan(sum(currPos,2)) &  (Dist2Centre < SET.BorderOutlierThreshold);
        OK.odor = currTheta>315 | currTheta<=45;
        OK.control = currTheta>135 & currTheta<=225;
        OK.Outside1 = currTheta>45 & currTheta<=135;
        OK.Outside2 = currTheta>225 & currTheta<=315;
        % Calculate the (proportion) of presence
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Presence.Quadrant.Odor(iTrial,1) = sum(OK.general & OK.odor);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Presence.Quadrant.Control(iTrial,1) = sum(OK.general & OK.control);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Presence.Quadrant.Outside1(iTrial,1) = sum(OK.general & OK.Outside1);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Presence.Quadrant.Outside2(iTrial,1) = sum(OK.general & OK.Outside2);
        % ---
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Odor(iTrial,1) = sum(OK.general & OK.odor) / sum(OK.general & (OK.odor | OK.control | OK.Outside1 | OK.Outside2));
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Control(iTrial,1) = sum(OK.general & OK.control) / sum(OK.general & (OK.odor | OK.control | OK.Outside1 | OK.Outside2));
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Outside1(iTrial,1) = sum(OK.general & OK.Outside1) / sum(OK.general & (OK.odor | OK.control | OK.Outside1 | OK.Outside2));
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Outside2(iTrial,1) = sum(OK.general & OK.Outside2) / sum(OK.general & (OK.odor | OK.control | OK.Outside1 | OK.Outside2));
        % Calculate the preference index
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).QuarterIndex(iTrial,1) = ...
            (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Odor(iTrial,1) - PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Control(iTrial,1)) / ...
            (PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Odor(iTrial,1) + PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Control(iTrial,1));
        
        
        
        % --- binnedPropQuarterTime_Odor
        clear OK
        for iAni = 1:DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).N
            currPos = [DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Pos_X(:,iAni), DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Pos_Y(:,iAni)];
            currTheta = atan2d(currPos(:,2), currPos(:,1));
            currTheta(currTheta<0) = 360-abs(currTheta(currTheta<0));
            % Get position to centre
            Dist2Centre = sqrt(sum(currPos'.*currPos'))';
            OK.odor(:,iAni) = currTheta>315 | currTheta<=45;
            OK.control(:,iAni) = currTheta>135 & currTheta<=225;
            OK.general(:,iAni) = ~isnan(sum(currPos,2)) &  (Dist2Centre < SET.BorderOutlierThreshold);
            OdorData(:,iAni) = (OK.general(:,iAni) & OK.odor(:,iAni));
            ControlData(:,iAni) = (OK.general(:,iAni) & OK.control(:,iAni));
            AnyData(:,iAni) = (OK.general(:,iAni) & (OK.odor(:,iAni) | OK.control(:,iAni)));
        end
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        if size(OdorData,2)>1
            OdorData = sum(OdorData,2)/size(OdorData,2);
        end
        OdorData = cumsum(OdorData)/size(OdorData,1);
        % Bin data into X bins
        [~,idx] = histc(vec, SET.BinVec);
        OdorData_bin = accumarray(idx(:), OdorData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Odor = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Odor; ...
            OdorData_bin'];
        
        % --- binnedPropQuarterTime_Control
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        if size(ControlData,2)>1
            ControlData = sum(ControlData,2)/size(ControlData,2);
        end
        ControlData = cumsum(ControlData)/size(ControlData,1);
        % Bin data into X bins
        [~,idx] = histc(vec, SET.BinVec);
        ControlData_bin = accumarray(idx(:), ControlData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Control = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Control; ...
            ControlData_bin'];
        
        % --- binnedPropQuarterTime_Any
        % Get time vector
        vec = 1:SET.CutAfter;
        % Get data
        if size(AnyData,2)>1
            AnyData = sum(AnyData,2)/size(AnyData,2);
        end
        AnyData = cumsum(AnyData)/size(AnyData,1);
        % Bin data into X bins
        [~,idx] = histc(vec, SET.BinVec);
        AnyData_bin = accumarray(idx(:), AnyData, [], @mean);
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Any = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Any; ...
            AnyData_bin'];
        
        % --- binnedPropQuarterTime_Delta
        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Delta = [...
            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Delta; ...
            OdorData_bin'-ControlData_bin'];
        clear vec OdorData* ControlData* idx OK
        
        clear Quadrant currTheta currPos Dist2Centre OK helper temp
    end%iTrial
    
    
    
    
    % ******************** Save averages in .csv ********************
    
    for iTrial = 1:size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Dist2Both,1)
        Pos_Inner = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).(['Trial_', num2str(iTrial)]).Pos_all;
        idx_Inner = (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > 1);
        Pos_Inner(idx_Inner,:) = [];
        avgVec = nanmedian([Pos_Inner(:,1), Pos_Inner(:,2)]);
        avgVec = avgVec/norm(avgVec);
        
        AvgData.ID{counterTable,1} = DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ID{iTrial,1};
        AvgData.Condition{counterTable,1} = SET.PlorOrder{cntCond}(1:3);
        AvgData.Odor_Concentration{counterTable,1} = SET.PlorOrder{cntCond}(5:end);
        AvgData.Condition_Odor_Concentration{counterTable,1} = SET.PlorOrder{cntCond};
        AvgData.PolarLocation_x(counterTable,1) = avgVec(1);
        AvgData.PolarLocation_y(counterTable,1) = avgVec(2);
        AvgData.StayDur_Odor(counterTable,1) =                      PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor(iTrial);
        AvgData.StayDur_Control(counterTable,1) =                   PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control(iTrial);
        AvgData.PropTimeSheltered_Odor(counterTable,1) =            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor(iTrial);
        AvgData.PropTimeSheltered_Control(counterTable,1) =         PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control(iTrial);
        AvgData.ShelterTimeIndex(counterTable,1) =                  PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex(iTrial);
        AvgData.EntryShelter_Odor(counterTable,1) =                 PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Odor(iTrial);
        AvgData.EntryShelter_Control(counterTable,1) =              PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryShelter_Control(iTrial);
        AvgData.EntryIndex(counterTable,1) =                        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex(iTrial);
        AvgData.ExitIndex(counterTable,1) =                         PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitIndex(iTrial);
        AvgData.PropEntryShelter_Odor(counterTable,1) =             PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Odor(iTrial);
        AvgData.PropEntryShelter_Control(counterTable,1) =          PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Control(iTrial);
        AvgData.ExitShelter_Odor(counterTable,1) =                  PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Odor(iTrial);
        AvgData.ExitShelter_Control(counterTable,1) =               PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ExitShelter_Control(iTrial);
        AvgData.PropExitShelter_Odor(counterTable,1) =              PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Odor(iTrial);
        AvgData.PropExitShelter_Control(counterTable,1) =           PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropExitShelter_Control(iTrial);
        AvgData.cumPropTimeSheltered_Odor(counterTable,1) =         PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Odor(iTrial,end);
        AvgData.cumPropTimeSheltered_Control(counterTable,1) =      PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Control(iTrial,end);
        AvgData.cumPropTimeSheltered_Delta(counterTable,1) =        PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).cumPropTimeSheltered_Delta(iTrial,end);
        AvgData.normcumPropTimeSheltered_Odor(counterTable,1) =     PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Odor(iTrial,end);
        AvgData.normcumPropTimeSheltered_Control(counterTable,1) =  PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Control(iTrial,end);
        AvgData.normcumPropTimeSheltered_Delta(counterTable,1) =    PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).normcumPropTimeSheltered_Delta(iTrial,end);
        AvgData.binnedPropTimeSheltered_Delta(counterTable,:) =     nanmean(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta);
        AvgData.SurvivalThreshold_Odor(counterTable,1) =            PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Odor(iTrial);
        AvgData.SurvivalThreshold_Control(counterTable,1) =         PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Control(iTrial);
        AvgData.SurvivalIndex(counterTable,1) =                     PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex(iTrial);
        AvgData.QuarterIndex(counterTable,1) =                         PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).QuarterIndex(iTrial);
        clear Pos_Inner idx_Inner avgVec avgTheta
        
        counterTable = counterTable+1;
        
    end%iTrial
    
    
    % Get min and max count of heatmaps to later cut at the minimum of the
    % max count and the maximum of the min count
    tempPos = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
    currCount = reshape(hist3(tempPos, [SET.HeatMapGrid , SET.HeatMapGrid]) / size(tempPos,1), [size(hist3(tempPos, [SET.HeatMapGrid , SET.HeatMapGrid]) / size(tempPos,1),1) * size(hist3(tempPos, [SET.HeatMapGrid , SET.HeatMapGrid]) / size(tempPos,1),2),1]);
    Q1 = quantile(currCount, 0.25);
    Q3 = quantile(currCount,0.75);
    Spread = 3*(Q3-Q1);
    MaxValue = Q3 + Spread;
    MinValue = Q1 - Spread;
    if MinValue<0; MinValue = 0; end
    maxCount = [maxCount; MaxValue];
    minCount = [minCount; MinValue];
    clear tempPos currCount Q1 Q3 Spread MaxValue MinValue
    
    
end%cntCond

% Save everything
save('tempDATA.mat', 'DATA', 'PooledDATA', 'AvgData', 'ANNOTATION', 'SET', 'maxSurvival', 'maxCount', 'minCount', '-v7.3')


%%                       ----- PLOT DATA -----
%**************************************************************************


% Already create figure windows and step-wise fill them
if SET.DoPlot
    % Heatmaps will be created on the fly for each condition. There is no
    % need to create the figure windows here.
    hFig.Polar_GRP_IND =                figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Polar_GRP_IND');
    hFig.Polar_Van_FecVan =             figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Polar_Van_FecVan');
    hFig.Polar_all =                    figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Polar_all');
    hFig.Dist2Both =                    figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Dist2Both');
    hFig.PropEntry =                    figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'PropEntry');
    hFig.RetentionTime =                figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'RetentionTime');
    hFig.cumTimeSheltered =             figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'cumTimeSheltered');
    hFig.cumTimeQuartered =             figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'cumTimeQuartered');
    hFig.PropTimeSheltered =            figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'PropTimeSheltered');
    hFig.Survival =                     figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Survival');
    hFig.Dendrogram =                   figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Dendrogram');
    hFig.SocialAttraction =             figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'SocialAttraction');
    hFig.PropPresence =                 figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'PropPresence');
    hFig.PreferenceIndices =            figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'PreferenceIndices');
    hFig.Speed =                        figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1], 'Name', 'Speed');
    
    xTick = [];
    
    % Iterate over all conditions
    for cntPlot = 1:length(SET.PlorOrder)
        
        
        % Get the right combination of condition, odor, and concentration
        currPlot = SET.PlorOrder(cntPlot,:);
        idx = cell2mat(strfind(currPlot, ' '));
        iCondition = find(strcmp(SET.ConditionNames, currPlot{1}(1:idx(1)-1)));
        iOdor = find(strcmp(SET.OdorNames, currPlot{1}(idx(1)+1:idx(2)-1)));
        iConcentration = find(strcmp(SET.ConcentrationNames, currPlot{1}(idx(2)+1:end)));
        clear idx currPlot
        
        
        % Create a string to identify the current concentration
        Concentration = ['C', SET.ConcentrationNames{iConcentration}];
        Concentration(Concentration=='.') = '_';
        
        
        % Continue, if data of the current condition, odor, and concentration have been pooled previously
        if ~isfield(DATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}), Concentration)
            error('Combination of condition, odor, and concentration does not exist')
        end
        % ************************* PLOTS *************************
        if SET.DoPlot
            xTick = [xTick;
                cntPlot];
            xTickLabel{1,cntPlot} = [SET.ConditionNames{iCondition}, ' | ', SET.OdorNames{iOdor}, ' | ', SET.ConcentrationNames{iConcentration}];
            xTickLabel{1,cntPlot}(xTickLabel{1,cntPlot} == '_') = '-';
            
            % Heatmap (divided into two zones: an inner, and an outer (border effects))
            % -------------------------------------------------------------
            hHeat=figure;
            set(gcf, 'units', 'centimeters', 'innerposition', [3 3 5 5])
            
            % Devide by inner and outer data to differentiate
            % between wall-fallowing and shelter-related
            % behaviour
            % --- Inner
            Pos_Inner = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
            idx_Inner = (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > 1);
            Pos_Inner(idx_Inner,:) = [];
            Pos_Inner = [Pos_Inner;...
                -1 -1
                -1 +1
                +1 -1
                +1 +1];
            % --- Outer
            Pos_Outer = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
            idx_Outer = (sqrt(sum(Pos_Outer'.*Pos_Outer'))' <= SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Outer'.*Pos_Outer'))' > 1);
            Pos_Outer(idx_Outer,:) = [];
            Pos_Outer = [Pos_Outer;...
                -1 -1
                -1 +1
                +1 -1
                +1 +1];
            % --- All
            Pos_All = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
            idx_All = (sqrt(sum(Pos_All'.*Pos_All'))' > 1);
            Pos_All(idx_All,:) = [];
            Pos_All = [Pos_All;...
                -1 -1
                -1 +1
                +1 -1
                +1 +1];
            
            % 2D hist of positions
            % --- Inner
            cntPresence_Inner = hist3(Pos_Inner, [SET.HeatMapGrid , SET.HeatMapGrid])/numel(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all);
            %         cntPresence_Inner=log10(cntPresence_Inner);
            %         cntPresence_Inner(isinf(cntPresence_Inner)) = 0;
            % --- Outer
            cntPresence_Outer = hist3(Pos_Outer, [SET.HeatMapGrid , SET.HeatMapGrid])/numel(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all);
            %         cntPresence_Outer=log10(cntPresence_Outer);
            %         cntPresence_Outer(isinf(cntPresence_Outer)) = 0;
            % --- All
            cntPresence_All = hist3(Pos_All, [SET.HeatMapGrid , SET.HeatMapGrid])/numel(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all);
            %         cntPresence_Inner=log10(cntPresence_Inner);
            %         cntPresence_Inner(isinf(cntPresence_Inner)) = 0;
            
            
            % Cut-off outliers
            cntPresence_Inner(cntPresence_Inner>min(maxCount)) = min(maxCount);
            cntPresence_Inner(cntPresence_Inner<max(minCount)) = max(minCount);
            cntPresence_Outer(cntPresence_Outer>min(maxCount)) = min(maxCount);
            cntPresence_Outer(cntPresence_Outer<max(minCount)) = max(minCount);
            cntPresence_All(cntPresence_Outer>min(maxCount)) = min(maxCount);
            cntPresence_All(cntPresence_Outer<max(minCount)) = max(minCount);
            
            cntPresence_Outer(cntPresence_Outer==0) = NaN;
                        
            
            % Plot everything together
            s = surf(cntPresence_All);
            s.EdgeColor = 'none';
            set(gca, 'view', [90 90])
            set(gca, 'XTick', [], 'YTick', [])
            axis equal
            xlim([1 , SET.HeatMapGrid])
            ylim([1 , SET.HeatMapGrid])
            caxis([round(max(minCount),5), round(min(maxCount),5)])
            colormap(viridis(50))
            caxis([round(max(minCount),5), round(min(maxCount),5)])
            
            % save
            export_fig(['FIG\raw\Heatmap_', SET.ConditionNames{iCondition}, '_', SET.OdorNames{iOdor}, '_', SET.ConcentrationNames{iConcentration}], '-png', '-r300')
            
            
            close(hHeat)
            % -------------------------------------------------------------
            
            
            
            
            % Heatmap relative density plots
            % -------------------------------------------------------------
            
            % --- First, divided by inner and outer zone ---
            hHeat_aux = figure;
            set(gcf, 'units', 'normalized', 'position', [0 0 1 1])
            
            subplot(1,2,1)
            hold on
            [cntPresence_Inner,c_Inner] = hist(Pos_Inner(:,1), SET.HeatMapGrid);
            [cntPresence_Outer,c_Outer] = hist(Pos_Outer(:,1), SET.HeatMapGrid);
            [cntPresence_All,~] = hist(Pos_All(:,1), SET.HeatMapGrid);
            plot(c_Inner, cntPresence_Inner/max(cntPresence_All), 'Color', SET.Color(cntPlot,:))
            plot(c_Outer, cntPresence_Outer/max(cntPresence_All), ':', 'Color', SET.Color(cntPlot,:))
            set(gca, 'units', 'centimeters', 'position', [1 5 7 2])
            set(gca, 'XTick', [], 'YTick', [0 1])
            xlim([-1 1])
            ylim([0 1])
            box on
            
            subplot(1,2,2)
            hold on
            
            [cntPresence_Inner,c_Inner] = hist(Pos_Inner(:,2), SET.HeatMapGrid);
            [cntPresence_Outer,c_Outer] = hist(Pos_Outer(:,2), SET.HeatMapGrid);
            [cntPresence_All,~] = hist(Pos_All(:,2), SET.HeatMapGrid);
            plot(cntPresence_Inner/max(cntPresence_All), fliplr(c_Inner), 'Color', SET.Color(cntPlot,:))
            plot(cntPresence_Outer/max(cntPresence_All), fliplr(c_Outer), ':', 'Color', SET.Color(cntPlot,:))
            set(gca, 'units', 'centimeters', 'position', [10 5 2 7])
            set(gca, 'XTick', [0 1], 'YTick', [])
            ylim([-1 1])
            xlim([0 1])
            box on
            
            export_fig(['FIG\raw\Heatmap_', SET.ConditionNames{iCondition}, '_', SET.OdorNames{iOdor}, '_', SET.ConcentrationNames{iConcentration},'_aux'], '-pdf')
            close(hHeat_aux)
            
            % clean-up
            clear cntPresence* Pos_All Pos_Inner Pos_Outer ax1 ax2 cb1 cb2 s c hHeat* idx_*
            % -------------------------------------------------------------
            
            
            
            
            % Polar plots
            % -------------------------------------------------------------
            if strcmp(SET.ConditionNames{iCondition}, 'GRP') || (strcmp(SET.ConditionNames{iCondition}, 'IND') && strcmp(SET.OdorNames{iOdor}, 'VAN') && strcmp(Concentration, 'C0_01'))
                figure(hFig.Polar_GRP_IND)
                % Avg over trials. Therefore, prepecify variables
                
                Pos_Inner = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
                idx_Inner = (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > 1);
                Pos_Inner(idx_Inner,:) = [];
                helper = Pos_Inner;
                helper = helper./sqrt(sum(helper'.*helper'))';
                theta = -atan2(helper(:,2), helper(:,1));
                h = polarhistogram(theta,11);
                hold on
                % Hist
                h.Normalization = 'probability';
                h.EdgeColor = 'none';
                h.FaceColor = SET.Color(cntPlot,:);
                h.FaceAlpha = 1;
                % Avg
                avgVec = nanmedian([Pos_Inner(:,1), Pos_Inner(:,2)]);
                v = polarplot([0, -atan2(avgVec(2), avgVec(1))], [0, norm(avgVec)], 'Color', SET.Color(cntPlot,:));
                v.LineWidth = 2;
                clear h v avgVec helper theta
                rlim([0 0.75])
            end
            
            if (strcmp(SET.ConditionNames{iCondition}, 'IND') && strcmp(SET.OdorNames{iOdor}, 'VAN') && strcmp(Concentration, 'C0_01'))...
                    || (strcmp(SET.ConditionNames{iCondition}, 'IND') && strcmp(SET.OdorNames{iOdor}, 'FEC_VAN') && strcmp(Concentration, 'C0_01_0_01'))
                figure(hFig.Polar_Van_FecVan)
                % Avg over trials. Therefore, prepecify variables
                Pos_Inner = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
                idx_Inner = (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > 1);
                Pos_Inner(idx_Inner,:) = [];
                helper = Pos_Inner;
                helper = helper./sqrt(sum(helper'.*helper'))';
                theta = -atan2(helper(:,2), helper(:,1));
                h = polarhistogram(theta,11);
                hold on
                % Hist
                h.Normalization = 'probability';
                h.EdgeColor = 'none';
                h.FaceColor = SET.Color(cntPlot,:);
                h.FaceAlpha = 1;
                % Avg
                avgVec = nanmedian([Pos_Inner(:,1), Pos_Inner(:,2)]);
                v = polarplot([0, -atan2(avgVec(2), avgVec(1))], [0, norm(avgVec)], 'Color', SET.Color(cntPlot,:));
                v.LineWidth = 2;
                clear h v avgVec helper theta
                rlim([0 0.75])
                
            end
            
            figure(hFig.Polar_all)
            Pos_Inner = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Pos_all;
            idx_Inner = (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > SET.BorderOutlierThreshold) | (sqrt(sum(Pos_Inner'.*Pos_Inner'))' > 1);
            Pos_Inner(idx_Inner,:) = [];
            avgVec = nanmedian([Pos_Inner(:,1), Pos_Inner(:,2)]);
            v = polarplot([0, -atan2(avgVec(2), avgVec(1))], [0, norm(avgVec)], 'Color', SET.Color(cntPlot,:));
            v.LineWidth = 2;
            rlim([0 0.75])
            hold on
            clear avgVec v h theta helper Pos_Inner
            % -------------------------------------------------------------
            
            
            
            
            % Histogram of distance to shelter (-1: control | 1: odor)
            % -------------------------------------------------------------
            figure(hFig.Dist2Both)
            subplot(3,3,cntPlot)
            hold on
            binWidth = (SET.HistBins(2)-SET.HistBins(1));
            binHalfWidth = binWidth/2;
            for iBin = 1:length(SET.HistBins)
                properties.BoxFaceCol = SET.Color(cntPlot,:);
                properties.nanmedianCol = 'k';
                properties.nanmedianWidth = 1;
                properties.OutlierSymbol = '.';
                properties.OutlierSize = 8;
                properties.WhiskerWidth = 1;
                boxplot_advanced(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Dist2Both(:,iBin), SET.HistBins(iBin), binWidth*0.95, properties)
            end%iBin
            clear iBin cnt_Hist
            xlim([-1-binHalfWidth, 1+binHalfWidth])
            xlabel('Relative Distance to Shelter Center')
            ylabel('Relative Time')
            set(gca, 'XTick', [-1 -0.5 0 0.5 1], 'XTickLabel', {'Control', '-0.5', '0', '0.5', 'Odor'})
            title([SET.ConditionNames{iCondition}, ' | ', SET.OdorNames{iOdor}, ' | ', SET.ConcentrationNames{iConcentration}], 'Interpreter', 'none')
            set(gcf, 'Color', 'w')
            ylim([0 1])
            clear properties
            % -------------------------------------------------------------
            
            
            
            
            % Barplots depicting the proportion of entries (only odor since it's mirrored)
            % -------------------------------------------------------------
            figure(hFig.PropEntry)
            hold on;
            
            data = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropEntryShelter_Odor;
            rng(4243)
            [data_boot, boot_idx] = bootstrp(SET.N_Boot, @mean, data);

            % --- Barplots
            plot([cntPlot,cntPlot],[mean(data_boot), mean(data_boot)+std(data_boot)],'k')
            plot([cntPlot-0.25, cntPlot+0.25],[mean(data_boot)+std(data_boot), mean(data_boot)+std(data_boot)],'k')
            rectangle('Position', [cntPlot-0.25, 0, 0.5, mean(data_boot)], 'FaceColor', SET.Color(cntPlot,:), 'EdgeColor', 'none')
            ylabel('Proportion of entries')
            
            % Statistics
            % Get the data
            z = data;
            PredVal = 0.5;
            z_tilde = z - mean(z) + PredVal;
            % Sample from the joined sample. For this, use the same indices
            z_boot = z_tilde(boot_idx);
            % Apply test statistic
            stat_observed = SET.TestStat1(z, length(z), PredVal);
            stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
            % Compute p values
            STATS.PropEntry.pVal(cntPlot) = mean(stat_boot >= stat_observed);
            STATS.PropEntry.Comparison{1, cntPlot} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            
            clear data data_boot currnanmean currstd CI_upper CI_lower z z_* PredVal stat_*
            % -------------------------------------------------------------
            
            
            
            
            % Duration animals spent in the shelter
            % -------------------------------------------------------------
            figure(hFig.RetentionTime)
            hold on
            % Get original data
            data.Odor = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor;
            data.Control = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control;
            % Resample data
            rng(4243)
            data.Odor_boot = bootstrp(SET.N_Boot, @mean, data.Odor);
            rng(4243)
            data.Control_boot = bootstrp(SET.N_Boot, @mean, data.Control);
            
            % Get mean, std, and CI
            % --- Odor
            currnanmean.Odor_boot = nanmean(data.Odor_boot);
            currnanstd.Odor_boot = nanstd(data.Odor_boot);
            currnansem.Odor_boot = currnanstd.Odor_boot / sqrt(size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Odor,1));
            CI_upper.Odor_boot = quantile(data.Odor_boot, 0.975);
            CI_lower.Odor_boot = quantile(data.Odor_boot, 0.025);
            % --- Control
            currnanmean.Control_boot = nanmean(data.Control_boot);
            currnanstd.Control_boot = nanstd(data.Control_boot);
            currnansem.Control_boot = currnanstd.Control_boot / sqrt(size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).StayDur_Control,1));
            CI_upper.Control_boot = quantile(data.Control_boot, 0.975);
            CI_lower.Control_boot = quantile(data.Control_boot, 0.025);
            
            % Bar plot, Odor
            rectangle('Position', [cntPlot-0.25, 0, 0.2, currnanmean.Odor_boot], 'FaceColor', SET.Color(cntPlot,:), 'EdgeColor', 'none')
            plot([cntPlot-0.15 cntPlot-0.15], [currnanmean.Odor_boot - currnanstd.Odor_boot, currnanmean.Odor_boot + currnanstd.Odor_boot], 'Color', [0 0 0], 'LineWidth', 1)
            plot([cntPlot-0.15-0.075 cntPlot-0.15+0.075], [currnanmean.Odor_boot + currnanstd.Odor_boot, currnanmean.Odor_boot + currnanstd.Odor_boot], 'Color', [0 0 0], 'LineWidth', 1)
            plot([cntPlot-0.15-0.075 cntPlot-0.15+0.075], [currnanmean.Odor_boot - currnanstd.Odor_boot, currnanmean.Odor_boot - currnanstd.Odor_boot], 'Color', [0 0 0], 'LineWidth', 1)
            
            % Bar plot, Control
            rectangle('Position', [cntPlot+0.05, 0, 0.2, currnanmean.Control_boot], 'FaceColor', 'none', 'EdgeColor', SET.Color(cntPlot,:))
            plot([cntPlot+0.15 cntPlot+0.15], [currnanmean.Control_boot - currnanstd.Control_boot, currnanmean.Control_boot + currnanstd.Control_boot], 'Color', [0 0 0], 'LineWidth', 1)
            plot([cntPlot+0.15-0.075 cntPlot+0.15+0.075], [currnanmean.Control_boot + currnanstd.Control_boot, currnanmean.Control_boot + currnanstd.Control_boot], 'Color', [0 0 0], 'LineWidth', 1)
            plot([cntPlot+0.15-0.075 cntPlot+0.15+0.075], [currnanmean.Control_boot - currnanstd.Control_boot, currnanmean.Control_boot - currnanstd.Control_boot], 'Color', [0 0 0], 'LineWidth', 1)
            set(gcf, 'Color', 'w')
            ylabel('Retention time (s)')
            clear data data_boot currnanmean CI_upper CI_lower
            % -------------------------------------------------------------
            
            
            
            
            % Speed over time
            % -------------------------------------------------------------
            figure(hFig.Speed)
            subplot(1,2,1)
            title('Speed over time')
            hold on;
            xvec = linspace(0, (size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed,2)/SET.FrameRate)/60, size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed,2));
            plot(xvec, nanmedian(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed), 'Color', SET.Color(cntPlot,:),'LineWidth', 3)
            ylabel('Average speed (cm/s)')
            xlabel('Time (min)')
            xlim([0, 30])
            
            subplot(1,2,2)
            title('Speed over time (moving only)')
            hold on;
            xvec = linspace(0, (size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed,2)/SET.FrameRate)/60, size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed,2));
            plot(xvec, nanmedian(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Speed_movOnly), 'Color', SET.Color(cntPlot,:),'LineWidth', 3)
            ylabel('Average speed (cm/s)')
            xlabel('Time (min)')
            xlim([0, 30])
            
            set(gcf, 'Color', 'w')
            % -------------------------------------------------------------
            
            
            
            
            % Cumulative Shelter Time
            % -------------------------------------------------------------
            figure(hFig.cumTimeSheltered)
            set(gcf, 'Color', 'w')
            
            % Statistics
            % Get data
            rng(4243)
            bootAny = bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Any);
            length_any = size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Any,1);
            rng(4243)
            [bootDelta, idx_delta] = bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta);
            length_delta = size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Delta,1);
            z_all = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Odor;
            y_all = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropTime_Control;
            length_z = size(z_all,1);
            length_y = size(y_all,1);
            
            bootPVal = nan(1,SET.PropTimeBins);
            for iBin= 1:SET.PropTimeBins
                
                % Get sampled data
                z = z_all(:,iBin);
                y = y_all(:,iBin);
                % Get joined distributions to sample from
                x = [z;y];
                % Draw bootstrap samples with replacement
                rng(4243)
                x_boot = reshape(randsample(x, length(x)*SET.N_Boot, 1), [length(x),SET.N_Boot]);
                % Split resamples according to length of odor and control
                % samples
                z_boot = x_boot(1:length_z,:);
                y_boot = x_boot(length_z+1:end,:);
                % Evaluate the test statistics on the original and on the
                % bootstrap samples
                teststat_sample = SET.TestStat2(z,      y,      length_z, length_y);
                teststat_boot =   SET.TestStat2(z_boot, y_boot, length_z, length_y);
                % Approximate achieved significance level
                bootPVal(iBin) = mean(teststat_boot >= teststat_sample);
                
            end%iBin
            
            STATS.CumTimeSheltered.pVal(cntPlot,:) = bootPVal;
            STATS.CumTimeSheltered.Comparison{cntPlot,1} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            
            clear z y z_* y_* teststat_* length_* idx*
            
            % Plot
            subplot(1,5,1); hold on
            plot(mean(bootAny), 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            ylabel('Cumul. prop. shelter time')
            xlabel('Bin #')
            
            subplot(1,5,2); hold on
            plot([cntPlot, cntPlot],[mean(bootAny(:,end))+std(bootAny(:,end)), mean(bootAny(:,end))-std(bootAny(:,end))],'k')
            plot(cntPlot, mean(bootAny(:,end)), 'o','MarkerEdgeColor', 'none', 'MarkerFaceColor', SET.Color(cntPlot,:), 'Markersize', 5)
            
            clear properties
            
            % Plot
            subplot(1,5,3); hold on
            plot(mean(bootDelta), 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            for iBin= 1:SET.PropTimeBins
                if bootPVal(iBin)<0.05
                    currFaceColor = SET.Color(cntPlot,:);
                else
                    currFaceColor = 'w';
                end
                plot(iBin, mean(bootDelta(:,iBin)), 'o', 'MarkerEdgeColor', SET.Color(cntPlot,:), 'MarkerFaceColor', currFaceColor)
            end%iBin
            ylabel('Diff. in cumul. prop. shelter time')
            xlabel('Bin #')
            
            subplot(1,5,4); hold on
            plot([cntPlot, cntPlot],[mean(bootDelta(:,end))+std(bootDelta(:,end)), mean(bootDelta(:,end))-std(bootDelta(:,end))],'k')
            plot(cntPlot, mean(bootDelta(:,end)), 'o','MarkerEdgeColor', 'none', 'MarkerFaceColor', SET.Color(cntPlot,:), 'Markersize', 5)
            
            % Plot
            subplot(1,5,5); hold on
            plot(bootPVal, 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            ylabel('p-Value')
            xlabel('Bin #')
            % -------------------------------------------------------------
            
            
            
            
            % Probability of being in a quarter over time
            % -------------------------------------------------------------
            figure(hFig.cumTimeQuartered)
            set(gcf, 'Color', 'w')
            
            % Statistics
            % Get data
            rng(4243)
            bootAny = bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Any);
            length_any = size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Any,1);
            rng(4243)
            [bootDelta, idx_delta] = bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Delta);
            length_delta = size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Delta,1);
            z_all = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Odor;
            y_all = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).binnedPropQuarterTime_Control;
            length_z = size(z_all,1);
            length_y = size(y_all,1);
            
            bootPVal = nan(1,SET.PropTimeBins);
            for iBin= 1:SET.PropTimeBins
                
                % Get sampled data
                z = z_all(:,iBin);
                y = y_all(:,iBin);
                % Get joined distributions to sample from
                x = [z;y];
                % Draw bootstrap samples with replacement
                rng(4243)
                x_boot = reshape(randsample(x, length(x)*SET.N_Boot, 1), [length(x),SET.N_Boot]);
                % Split resamples according to length of odor and control
                % samples
                z_boot = x_boot(1:length_z,:);
                y_boot = x_boot(length_z+1:end,:);
                % Evaluate the test statistics on the original and on the
                % bootstrap samples
                teststat_sample = SET.TestStat2(z,      y,      length_z, length_y);
                teststat_boot =   SET.TestStat2(z_boot, y_boot, length_z, length_y);
                % Approximate achieved significance level
                bootPVal(iBin) = mean(teststat_boot >= teststat_sample);
                
            end%iBin
            
            STATS.cumTimeQuartered.pVal(cntPlot,:) = bootPVal;
            STATS.cumTimeQuartered.Comparison{cntPlot,1} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            
            clear z y z_* y_* teststat_* length_* idx*
            
            % Plot
            subplot(1,5,1); hold on
            plot(mean(bootAny), 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            ylabel('Cumul. prop. quarter time')
            xlabel('Bin #')
            
            subplot(1,5,2); hold on
            plot([cntPlot, cntPlot],[mean(bootAny(:,end))+std(bootAny(:,end)), mean(bootAny(:,end))-std(bootAny(:,end))],'k')
            plot(cntPlot, mean(bootAny(:,end)), 'o','MarkerEdgeColor', 'none', 'MarkerFaceColor', SET.Color(cntPlot,:), 'Markersize', 5)
            
            clear properties
            
            % Plot
            subplot(1,5,3); hold on
            plot(mean(bootDelta), 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            for iBin= 1:SET.PropTimeBins
                if bootPVal(iBin)<0.05
                    currFaceColor = SET.Color(cntPlot,:);
                else
                    currFaceColor = 'w';
                end
                plot(iBin, mean(bootDelta(:,iBin)), 'o', 'MarkerEdgeColor', SET.Color(cntPlot,:), 'MarkerFaceColor', currFaceColor)
            end%iBin
            ylabel('Diff. in cumul. prop. quarter time')
            xlabel('Bin #')
            
            subplot(1,5,4); hold on
            plot([cntPlot, cntPlot],[mean(bootDelta(:,end))+std(bootDelta(:,end)), mean(bootDelta(:,end))-std(bootDelta(:,end))],'k')
            plot(cntPlot, mean(bootDelta(:,end)), 'o','MarkerEdgeColor', 'none', 'MarkerFaceColor', SET.Color(cntPlot,:), 'Markersize', 5)
            
            % Plot
            subplot(1,5,5); hold on
            plot(bootPVal, 'Color', SET.Color(cntPlot,:), 'LineWidth', 2)
            ylabel('p-Value')
            xlabel('Bin #')
            % -------------------------------------------------------------
            
            
            
            
            figure(hFig.PropTimeSheltered)
            % -------------------------------------------------------------
            hold on
            % Get data
            % --- odor
            odorData = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Odor;
            rng(4243)
            [odorData_boot, idx_odor] = bootstrp(SET.N_Boot, @mean, odorData);
            odorData_avg = mean(odorData_boot);
            odorData_std = std(odorData_boot);
            % --- control
            controlData = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropTimeSheltered_Control;
            rng(4243)
            [controlData_boot, idx_control] = bootstrp(SET.N_Boot, @mean, controlData);
            controlData_avg = mean(controlData_boot);
            controlData_std = std(controlData_boot);
            
            
            % Plot means as bar plots
            rectangle('Position', [cntPlot*2, 0, 0.5, odorData_avg], 'EdgeColor', 'none', 'FaceColor', SET.Color(cntPlot,:))
            rectangle('Position', [cntPlot*2-0.5, 0, 0.5, controlData_avg], 'FaceColor', 'w', 'EdgeColor', SET.Color(cntPlot,:))
            % Indicate variability
            plot([cntPlot*2-0.25, cntPlot*2-0.25],[controlData_avg, controlData_avg + controlData_std], 'k')
            plot([cntPlot*2-0.35, cntPlot*2-0.15],[controlData_avg + controlData_std, controlData_avg + controlData_std], 'k')
            plot([cntPlot*2+0.25, cntPlot*2+0.25],[odorData_avg, odorData_avg + odorData_std], 'k')
            plot([cntPlot*2+0.15, cntPlot*2+0.35],[odorData_avg + odorData_std, odorData_avg + odorData_std], 'k')

            
            % Statistics
            % Get sampled data
            % --- odor
            z = odorData;
            length_z = length(z);
            % --- control
            y = controlData;
            length_y = length(y);
            % Get joined distribution to sample from
            x = [z;y];
            % Draw bootstrap samples with replacement
            rng(4243)
            x_boot = reshape(randsample(x, length(x)*SET.N_Boot, 1), [length(x),SET.N_Boot]);
            % Split resamples according to length of odor and control
            % samples
            z_boot = x_boot(1:length_z,:);
            y_boot = x_boot(length_z+1:end,:);
            % Evaluate the test statistics on the original and on the
            % bootstrap samples
            teststat_sample = SET.TestStat2(z,      y,      length_z, length_y);
            teststat_boot =   SET.TestStat2(z_boot, y_boot, length_z, length_y);
            % Approximate achieved significance level
            STATS.PropTimeSheltered.pVal(cntPlot) = mean(teststat_boot >= teststat_sample);
            STATS.PropTimeSheltered.Comparison{1, cntPlot} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            clear z y length_* z_* y_* teststat_* idx*
            clear odorData* controlData* out1Data* out2Data* t1 t2
            % -------------------------------------------------------------

            
            
            
            % Survival Curves
            % -------------------------------------------------------------
            xvec = linspace(0, (SET.CutAfter/SET.FrameRate)/60, SET.CutAfter);
            
            % Resample data at threshold
            SurvivalData_odor = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Odor;
            SurvivalData_odor(isnan(SurvivalData_odor)) = [];
            rng(4243)
            [BootSurvival_odor, idx_odor] = bootstrp(SET.N_Boot, @mean, SurvivalData_odor);
            
            SurvivalData_control = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalThreshold_Control;
            SurvivalData_control(isnan(SurvivalData_control)) = [];
            rng(4243)
            [BootSurvival_control, idx_control] = bootstrp(SET.N_Boot, @mean, SurvivalData_control);
            
            
            % --- Save individual survival figure
            tempFig = figure('units', 'normalized', 'position', [0 0 1 1], 'Color', 'w');
            hold on
            plot(xvec, nanmean(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor), 'Color', SET.Color(cntPlot,:), 'LineWidth', 1)
            plot(xvec, nanmean(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control),'Color', [0.85 0.85 0.85],  'LineWidth', 1)
            set(gca, 'units', 'centimeters', 'position', [5 5 0.75 1.5])
            xlim([xvec(1), xvec(end)])
            ylim([0 1])
            export_fig(['FIG\raw\survival_',SET.ConditionNames{iCondition},'_',SET.OdorNames{iOdor},'_',Concentration],'-pdf')
            close(tempFig); clear tempFig
            
            % --- Open figure
            figure(hFig.Survival)
            set(gcf, 'Color', 'w')
            
            % --- Survival in control shelter
            subplot(1,3,1)
            hold on
            plot(xvec, nanmean(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor), 'Color', SET.Color(cntPlot,:), 'LineWidth', 3)
            ylabel('Leaving probability')
            xlabel('Time (min)')
            xlim([0 2.5])
            set(gcf, 'Color', 'w')
            
            % --- Survival in odor shelter
            subplot(1,3,2)
            hold on
            xvec = linspace(0, (size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor,2)/SET.FrameRate)/60, size(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Odor,2));
            plot(xvec, nanmean(PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).Survival_Control), 'Color', SET.Color(cntPlot,:), 'LineWidth', 3)
            ylabel('Leaving probability')
            xlabel('Time (min)')
            xlim([0 2.5])
            
            % --- Difference in survival time at threshold
            subplot(1,3,3)
            hold on
            properties.NumPoints = 250;
            properties.EdgeCol = SET.Color(cntPlot,:);
            properties.OutlierCol = SET.Color(cntPlot,:);
            properties.OutlierSize = 1;
            properties.AvgType = 'mean';
            properties.MeanWidth = 1;
            properties.MinVal = 0;
            properties.MeanCol = 'k';
            violinplot_advanced(BootSurvival_odor, cntPlot-0.2, 0.15, properties)
            properties.MeanCol = SET.Color(cntPlot,:);
            violinplot_advanced(BootSurvival_control, cntPlot+0.2, 0.15, properties)
            clear properties
            xlim([0.5 8.5])
            
            % Statistics
            % Get sampled data
            z = SurvivalData_odor;
            length_z = length(z);
            y = SurvivalData_control;
            length_y = length(y);
            % Get joined distribution to sample from
            x = [z;y];
            % Draw bootstrap samples with replacement
            rng(4243)
            x_boot = reshape(randsample(x, length(x)*SET.N_Boot, 1), [length(x),SET.N_Boot]);
            % Split resamples according to length of odor and control
            % samples
            z_boot = x_boot(1:length_z,:);
            y_boot = x_boot(length_z+1:end,:);
            % Evaluate the test statistics on the original and on the
            % bootstrap samples
            teststat_sample = SET.TestStat2(z,      y,      length_z, length_y);
            teststat_boot =   SET.TestStat2(z_boot, y_boot, length_z, length_y);
            % Approximate achieved significance level
            STATS.Survival.pVal(cntPlot) = mean(teststat_boot >= teststat_sample);
            STATS.Survival.Comparison{1, cntPlot} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            clear z y length_* z_* y_* teststat_* idx* BootSurvival*
            % -------------------------------------------------------------
            
                        
            
            
            % PropPresence
            % -------------------------------------------------------------
            figure(hFig.PropPresence)
            hold on
            % Get data
            % --- odor
            odorData = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Odor;
            rng(4243)
            [odorData_boot, idx_odor] = bootstrp(SET.N_Boot, @mean, odorData);
            odorData_avg = mean(odorData_boot);
            odorData_std = std(odorData_boot);
            % --- control
            controlData = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Control;
            rng(4243)
            [controlData_boot, idx_control] = bootstrp(SET.N_Boot, @mean, controlData);
            controlData_avg = mean(controlData_boot);
            controlData_std = std(controlData_boot);
            % --- outside1
            out1Data = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Outside1;
            rng(4243)
            out1Data_boot = bootstrp(SET.N_Boot, @mean, out1Data);
            out1Data_avg = mean(out1Data_boot);
            out1Data_std = std(out1Data_boot);
            % --- outside2
            out2Data = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).PropPresence.Quadrant.Outside2;
            rng(4243)
            out2Data_boot = bootstrp(SET.N_Boot, @mean, out2Data);
            out2Data_avg = mean(out2Data_boot);
            out2Data_std = std(out2Data_boot);
            
            % Plot means as bar plots
            rectangle('Position', [cntPlot*2, 0, 0.5, odorData_avg], 'EdgeColor', 'none', 'FaceColor', SET.Color(cntPlot,:))
            rectangle('Position', [cntPlot*2-0.5, 0, 0.5, controlData_avg], 'FaceColor', 'w', 'EdgeColor', SET.Color(cntPlot,:))
            rectangle('Position', [cntPlot*2-0.5, -out1Data_avg, 0.5, out1Data_avg], 'EdgeColor', 'k', 'FaceColor', [0.85 0.85 0.85])
            rectangle('Position', [cntPlot*2, -out2Data_avg, 0.5, out2Data_avg], 'EdgeColor', 'k', 'FaceColor', [0.85 0.85 0.85])
            % Indicate variability
            plot([cntPlot*2-0.25, cntPlot*2-0.25],[controlData_avg, controlData_avg + controlData_std], 'k')
            plot([cntPlot*2-0.35, cntPlot*2-0.15],[controlData_avg + controlData_std, controlData_avg + controlData_std], 'k')
            plot([cntPlot*2+0.25, cntPlot*2+0.25],[odorData_avg, odorData_avg + odorData_std], 'k')
            plot([cntPlot*2+0.15, cntPlot*2+0.35],[odorData_avg + odorData_std, odorData_avg + odorData_std], 'k')
            plot([cntPlot*2-0.25, cntPlot*2-0.25],[-out1Data_avg, -out1Data_avg - out1Data_std], 'k')
            plot([cntPlot*2-0.35, cntPlot*2-0.15],[-out1Data_avg - out1Data_std, -out1Data_avg - out1Data_std], 'k')
            plot([cntPlot*2+0.25, cntPlot*2+0.25],[-out2Data_avg, -out2Data_avg - out2Data_std], 'k')
            plot([cntPlot*2+0.15, cntPlot*2+0.35],[-out2Data_avg - out2Data_std, -out2Data_avg - out2Data_std], 'k')
            
            
            % Statistics
            % Get sampled data
            % --- odor
            z = odorData;
            length_z = length(z);
            % --- control
            y = controlData;
            length_y = length(y);
            % Get joined distribution to sample from
            x = [z;y];
            % Draw bootstrap samples with replacement
            rng(4243)
            x_boot = reshape(randsample(x, length(x)*SET.N_Boot, 1), [length(x),SET.N_Boot]);
            % Split resamples according to length of odor and control
            % samples
            z_boot = x_boot(1:length_z,:);
            y_boot = x_boot(length_z+1:end,:);
            % Evaluate the test statistics on the original and on the
            % bootstrap samples
            teststat_sample = SET.TestStat2(z,      y,      length_z, length_y);
            teststat_boot =   SET.TestStat2(z_boot, y_boot, length_z, length_y);
            % Approximate achieved significance level
            STATS.PropPresence.pVal(cntPlot) = mean(teststat_boot >= teststat_sample);
            STATS.PropPresence.Comparison{1, cntPlot} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
            clear z y length_* z_* y_* teststat_* idx*
            clear odorData* controlData* out1Data* out2Data* t1 t2
            % -------------------------------------------------------------
            
            
            legend_str{1, cntPlot} = [SET.ConditionNames{iCondition}, ' | ', SET.OdorNames{iOdor}, ' | ', SET.ConcentrationNames{iConcentration}];
        end
    end%cntPlot
    
end


%% Save various variables to table. --> Could be used for GLMMs

T = table(...
    AvgData.ID,...
    AvgData.Condition,...
    AvgData.Odor_Concentration,...
    AvgData.PolarLocation_x,...
    AvgData.PolarLocation_y,...
    AvgData.StayDur_Odor,...
    AvgData.StayDur_Control,...
    AvgData.PropTimeSheltered_Odor,...
    AvgData.PropTimeSheltered_Control,...
    AvgData.ShelterTimeIndex,...
    AvgData.EntryShelter_Odor,...
    AvgData.EntryShelter_Control,...
    AvgData.PropEntryShelter_Odor,...
    AvgData.PropEntryShelter_Control,...
    AvgData.ExitShelter_Odor,...
    AvgData.ExitShelter_Control,...
    AvgData.PropExitShelter_Odor,...
    AvgData.PropExitShelter_Control,...
    AvgData.cumPropTimeSheltered_Control,...
    AvgData.cumPropTimeSheltered_Odor,...
    AvgData.normcumPropTimeSheltered_Control,...
    AvgData.normcumPropTimeSheltered_Odor,...
    AvgData.SurvivalThreshold_Control,...
    AvgData.SurvivalThreshold_Odor,...
    AvgData.SurvivalIndex,...
    AvgData.QuarterIndex,...
    'VariableNames',...
    {...
    'ID';...
    'Condition';...
    'Odor_Concentration';...
    'PolarLocation_x';...
    'PolarLocation_y';...
    'StayDur_Odor';...
    'StayDur_Control';...
    'PropTimeSheltered_Odor';...
    'PropTimeSheltered_Control';...
    'ShelterTimeIndex';...
    'EntryShelter_Odor';...
    'EntryShelter_Control';...
    'PropEntryShelter_Odor';...
    'PropEntryShelter_Control';...
    'ExitShelter_Odor';...
    'ExitShelter_Control';...
    'PropExitShelter_Odor';...
    'PropExitShelter_Control';...
    'cumPropTimeSheltered_Control';...
    'cumPropTimeSheltered_Odor';...
    'normcumPropTimeSheltered_Control';...
    'normcumPropTimeSheltered_Odor';...
    'SurvivalThreshold_Control';...
    'SurvivalThreshold_Odor';...
    'SurvivalIndex';...
    'QuarterIndex';....
    });
writetable(T,'AvgData.csv')
clc



%% Preference Indices

% Iterate over all conditions
figure(hFig.PreferenceIndices)
width = 0.5;
clear boot_idx
for cntCond = 1:length(SET.PlorOrder)
    
    % Get the right combination of condition, odor, and concentration
    currCond = SET.PlorOrder(cntCond,:);
    idx = cell2mat(strfind(currCond, ' '));
    iCondition = find(strcmp(SET.ConditionNames, currCond{1}(1:idx(1)-1)));
    iOdor = find(strcmp(SET.OdorNames, currCond{1}(idx(1)+1:idx(2)-1)));
    iConcentration = find(strcmp(SET.ConcentrationNames, currCond{1}(idx(2)+1:end)));
    clear idx currPlot    
    % Create a string to identify the current concentration
    Concentration = ['C', SET.ConcentrationNames{iConcentration}];
    Concentration(Concentration=='.') = '_';
    
    
    rng(4243)
    [ShelterTimeIndex_boot, boot_idx.ShelterTimeIndex] = bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex);
    rng(4243)
    [EntryIndex_boot,       boot_idx.EntryIndex] =       bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex);
    rng(4243)
    [QuarterIndex_boot,     boot_idx.QuarterIndex] =     bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).QuarterIndex);
    rng(4243)
    [SurvivalIndex_boot,    boot_idx.SurvivalIndex] =    bootstrp(SET.N_Boot, @mean, PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex);
    
    % --- ShelterTimeIndex_boot
    subplot(1,4,1); hold on
    properties.NumPoints = 250;
    properties.EdgeCol = SET.Color(cntCond,:);
    properties.OutlierCol = SET.Color(cntCond,:);
    properties.OutlierSize = 1;
    properties.MaxVal = 1;
    properties.MinVal = -1;
    properties.AvgType = 'mean';
    properties.MeanWidth = 1;
    violinplot_advanced(ShelterTimeIndex_boot, cntCond, width, properties)
    xlim([0.5 8.5]); ylim([-1 1])
    ylabel('ShelterTimeIndex')
    
    % Statistics vs group level
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex;
    PredVal = mean(PooledDATA.GRP.VAN.C0_01.ShelterTimeIndex);
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.ShelterTimeIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.ShelterTimeIndex_vsGRP.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.ShelterTimeIndex_vsGRP.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal    
    
    % Statistics vs zero
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).ShelterTimeIndex;
    PredVal = 0;
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.ShelterTimeIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.ShelterTimeIndex_vs0.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.ShelterTimeIndex_vs0.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
    
    
    
    
    % --- EntryIndex_boot
    subplot(1,4,2); hold on
    properties.NumPoints = 250;
    properties.EdgeCol = SET.Color(cntCond,:);
    properties.OutlierCol = SET.Color(cntCond,:);
    properties.OutlierSize = 1;
    properties.MaxVal = 1;
    properties.MinVal = -1;
    properties.AvgType = 'mean';
    properties.MeanWidth = 1;
    violinplot_advanced(EntryIndex_boot, cntCond, width, properties)
    xlim([0.5 8.5]); ylim([-1 1])
    ylabel('EntryIndex')    
    
    % Statistics vs group level
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex;
    PredVal = mean(PooledDATA.GRP.VAN.C0_01.EntryIndex);
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.EntryIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.EntryIndex_vsGRP.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.EntryIndex_vsGRP.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
        
    % Statistics vs zero
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).EntryIndex;
    PredVal = 0;
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.EntryIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.EntryIndex_vs0.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.EntryIndex_vs0.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
    
    
    
    
    % --- QuarterIndex_boot
    subplot(1,4,3); hold on
    properties.NumPoints = 250;
    properties.EdgeCol = SET.Color(cntCond,:);
    properties.OutlierCol = SET.Color(cntCond,:);
    properties.OutlierSize = 1;
    properties.MaxVal = 1;
    properties.MinVal = -1;
    properties.AvgType = 'mean';
    properties.MeanWidth = 1;
    violinplot_advanced(QuarterIndex_boot, cntCond, width, properties)
    xlim([0.5 8.5]); ylim([-1 1])
    ylabel('QuarterIndex')
    
    % Statistics vs group level
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).QuarterIndex;
    PredVal = mean(PooledDATA.GRP.VAN.C0_01.QuarterIndex);
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.QuarterIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.QuarterIndex_vsGRP.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.QuarterIndex_vsGRP.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
    
    % Statistics  vs zero
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).QuarterIndex;
    PredVal = 0;
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.QuarterIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.QuarterIndex_vs0.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.QuarterIndex_vs0.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
    
    
    
    
    % --- SurvivalIndex_boot
    subplot(1,4,4); hold on
    properties.NumPoints = 250;
    properties.EdgeCol = SET.Color(cntCond,:);
    properties.OutlierCol = SET.Color(cntCond,:);
    properties.OutlierSize = 1;
    properties.MaxVal = 1;
    properties.MinVal = -1;
    properties.AvgType = 'mean';
    properties.MeanWidth = 1;
    violinplot_advanced(SurvivalIndex_boot, cntCond, width, properties)
    xlim([0.5 8.5]); ylim([-1 1])
    ylabel('SurvivalIndex')
    
    % Statistics vs group level
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex;
    PredVal = mean(PooledDATA.GRP.VAN.C0_01.SurvivalIndex);
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.SurvivalIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.SurvivalIndex_vsGRP.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.SurvivalIndex_vsGRP.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
    
    % Statistics vs group level
    % Get the data
    z = PooledDATA.(SET.ConditionNames{iCondition}).(SET.OdorNames{iOdor}).(Concentration).SurvivalIndex;
    PredVal = 0;
    z_tilde = z - mean(z) + PredVal;
    % Sample from the joined sample. For this, use the same indices
    z_boot = z_tilde(boot_idx.SurvivalIndex);
    % Apply test statistic
    stat_observed = SET.TestStat1(z, length(z), PredVal);
    stat_boot = SET.TestStat1(z_boot, length(z), PredVal);
    % Compute p values
    STATS.SurvivalIndex_vs0.pVal(cntCond) = mean(stat_boot >= stat_observed);
    STATS.SurvivalIndex_vs0.Comparison{1, cntCond} = [SET.ConditionNames{iCondition}, '-',SET.OdorNames{iOdor}, '-',Concentration];
    clear z z_tilde z_boot stat_observed stat_boot PredVal
        
    clear properties
end






%% Dendrogram based on behavioural data


uniqueNames = unique(AvgData.Condition_Odor_Concentration);
clear X X_boot
for iCond = 1:length(uniqueNames)
    
    idx = find(strcmp(AvgData.Condition_Odor_Concentration, uniqueNames(iCond)));
    
    X(iCond,:) = [...
        nanmean(AvgData.ShelterTimeIndex(idx)),...
        nanmean(AvgData.QuarterIndex(idx)),...
        nanmean(AvgData.EntryIndex(idx)),...
        nanmean(AvgData.SurvivalIndex(idx)),...
        ];

end

figure(hFig.Dendrogram)

Y = pdist(X,'euclidean');
Z = linkage(Y,'average');
leafOrder = optimalleaforder(Z,Y);
dendrogram(Z, 'Reorder', leafOrder)
% set(gca, 'YTick', 0:7)
set(gca, 'XTickLabels', uniqueNames(str2num(get(gca, 'XTickLabel'))), 'view', [90 90])
STATS.Dendrogram = cophenet(Z,Y);
hFig.Dendrogram = gcf;


%% Social Attraction
figure(hFig.SocialAttraction)
LocNames = {'Odor','Control','Outside1','Outside2'};
LocCol = [...
    112 174 224;...
    216.75 216.75 216.75;...
    0 0 0;...
    0 0 0]/255;

% Inter-animal attraction I
% -------------------------------------------------------------------------
% (Time bout duration VS number of animals at the respective location
% during stay)

cnt_Odor=1;
cnt_Control=1;
clear corrOdor corrControl
PopLoc.Odor_all = [];
PopLoc.Control_all = [];
for iTrial = 1:length(fields(DATA.GRP.VAN.C0_01))-4
    
    % Get how many animals were sheltered at point in time
    PopLoc.Odor = sum(DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Odor');
    PopLoc.Control = sum(DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Control');
    
    % Combine everything for histogram
    PopLoc.Odor_all = [PopLoc.Odor_all, PopLoc.Odor];
    PopLoc.Control_all = [PopLoc.Control_all, PopLoc.Control];
    
    for iAni = 1:size(DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Odor, 2)
        
        % Get the shelter time for each animal
        AniLoc.Odor = DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Odor(:,iAni);
        AniLoc.Control = DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Control(:,iAni);
        
        % Label stays
        Odor = bwlabel(AniLoc.Odor);
        Control = bwlabel(AniLoc.Control);
        
        % Get list of stays
        unique_Odor = unique(Odor);
        unique_Control = unique(Control);
        
        % Iterate over stays in odor shelter
        for iClust = (unique_Odor(unique_Odor~=0))'
            corrOdor(1,cnt_Odor) = round(nanmean(PopLoc.Odor(find(Odor==iClust))),0);
            corrOdor(2,cnt_Odor) = sum(Odor==iClust);
            cnt_Odor = cnt_Odor+1;
        end%iClust
        
        % Iterate over stays in control shelter
        for iClust = (unique_Control(unique_Control~=0 & unique_Control~=unique_Control(end)))'
            corrControl(1,cnt_Control) = round(nanmean(PopLoc.Control(find(Control==iClust))),0);
            corrControl(2,cnt_Control) = sum(Control==iClust);
            cnt_Control = cnt_Control+1;
        end%iClust
        
    end%iAni
end%iTrial


% Histogram
% ---------
subplot(1,3,1)
xVec_Odor = min(PopLoc.Odor_all):max(PopLoc.Odor_all);
xVec_Control = min(PopLoc.Control_all):max(PopLoc.Control_all);
for iBin = 1:length(xVec_Odor)
    c1(iBin) = sum((PopLoc.Odor_all == xVec_Odor(iBin))) / length(PopLoc.Odor_all);
end
for iBin = 1:length(xVec_Control)
    c2(iBin) = sum((PopLoc.Control_all == xVec_Control(iBin))) / length(PopLoc.Control_all);
end


hold on
% Hist for odor shelter
for iBin = 1:length(c1)
    rectangle('Position', [xVec_Odor(iBin)-(0.5*0.75), 0, 0.5, c1(iBin)*100], 'FaceColor', LocCol(find(strcmp(LocNames, 'Odor')),:), 'EdgeColor', 'none')
end
% Hist for control shelter
for iBin = 1:length(c2)
    rectangle('Position', [xVec_Control(iBin)-(0.5*0.75)+0.25, 0, 0.5, c2(iBin)*100], 'FaceColor', LocCol(find(strcmp(LocNames, 'Control')),:), 'EdgeColor', 'none')
end
xlim([min([xVec_Odor, xVec_Control])-1, max([xVec_Odor, xVec_Control])+1])
set(gca, 'XTick', unique([xVec_Odor, xVec_Control]))

% Calculate a chi-squared test, comparing odor/control & alone/not alone
% --- Observed data
n1 = round(c1(2)*100, 0); N1 = round(sum(c1(3:end)*100), 0);
n2 = round(c2(2)*100, 0); N2 = round(sum(c2(3:end)*100), 0);
% --- Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2);
% --- Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% --- Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected);
p = 1 - chi2cdf(chi2stat,1);
% --- Save results
STATS.SocialHistogram.pVal = p;
STATS.SocialHistogram.ChiVal = chi2stat;
STATS.SocialHistogram.ChiComp = [n1,N1; n2,N2];


% Inter-animal attraction I
subplot(1,3,2)
hold on
cnt=1;
iAniVector = round(unique(corrOdor(1,:)),0);
for iAni = iAniVector
    
    if iAni >= 6
        idx = find(corrOdor(1,:) > iAniVector(end-3));
        QuitAfter = 1;
    else
        QuitAfter = 0;
        idx = find(corrOdor(1,:) == iAni);
    end
    currData = log10(corrOdor(2,idx)/SET.FrameRate);
    if length(currData)>1
        rng(4243)
        boot_data = bootstrp(SET.N_Boot, @mean, currData);
        boot_mean = mean(boot_data);
        boot_std = std(boot_data);
    else
        boot_data = currData;
        boot_mean = currData;
        boot_std = 0;
    end
    
    %--- Beeswarmplot
    %     properties.MarkerFaceColor = LocCol(find(strcmp(LocNames, 'Odor')),:);
    %     properties.MarkerSize = 2;
    %     beeswarmplot_advanced(boot_data, iAni, 0.5, properties)
    %     clear properties
    %     plot([iAni-0.25, iAni+0.25],[mean(currData), mean(currData)], 'Color', LocCol(find(strcmp(LocNames, 'Odor')),:), 'LineWidth', 2)
    
    properties.NumPoints = 250;
    properties.EdgeCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    properties.MeanWidth = 1;
    properties.OutlierCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    properties.OutlierSize = 1;
    properties.AvgType = 'mean';
    properties.MeanCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    violinplot_advanced(boot_data, iAni, 0.5, properties)
    clear properties
    
    %     IAA_1.odor.x(cnt,:) = ones(1,length(boot_data))*iAni;
    %     IAA_1.odor.y(cnt,:) = boot_data;
    %     IAA_1.odor.avg_y(cnt,1) = boot_mean;
    IAA_1.odor.N(cnt,1) = length(idx);
    cnt=cnt+1;
    
    if QuitAfter
        clear QuitAfter
        break
    end
end



cnt=1;
iAniVector = round(unique(corrControl(1,:)),0);
for iAni = iAniVector
    
    if iAni >= 6
        idx = find(corrControl(1,:) > iAniVector(end-3));
        QuitAfter = 1;
    else
        QuitAfter = 0;
        idx = find(corrControl(1,:) == iAni);
    end
    currData = log10(corrControl(2,idx)/SET.FrameRate);
    if length(currData)>1
        rng(4243)
        boot_data = bootstrp(SET.N_Boot, @mean, currData);
        boot_mean = mean(boot_data);
        boot_std = std(boot_data);
    else
        boot_data = currData;
        boot_mean = currData;
        boot_std = 0;
    end
    
    %--- Beeswarmplot
    %     properties.MarkerFaceColor = LocCol(find(strcmp(LocNames, 'Control')),:);
    %     properties.MarkerSize = 2;
    %     beeswarmplot_advanced(boot_data, iAni, 0.5, properties)
    %     clear properties
    %     plot([iAni-0.25, iAni+0.25],[mean(currData), mean(currData)], 'Color', LocCol(find(strcmp(LocNames, 'Control')),:), 'LineWidth', 2)
    
    properties.NumPoints = 250;
    properties.EdgeCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    properties.MeanWidth = 1;
    properties.OutlierCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    properties.OutlierSize = 1;
    properties.AvgType = 'mean';
    properties.MeanCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    violinplot_advanced(boot_data, iAni, 0.5, properties)
    clear properties
    
    %     IAA_1.control.x(cnt,:) = ones(1,length(boot_data))*iAni;
    %     IAA_1.control.y(cnt,:) = boot_data;
    %     IAA_1.control.avg_y(cnt,1) = boot_mean;
    IAA_1.control.N(cnt,1) = length(idx);
    cnt=cnt+1;
    
    if QuitAfter
        clear QuitAfter
        break
    end
    
end

% ylabel('Time bount in shelter (log10(s))')
% xlabel('Avg. # animals sheltered during stay')
set(gca, 'XTick', round(unique(corrControl(1,:)),0))
xlim([min(round(unique(corrControl(1,:)),0))-1, max(round(unique(corrControl(1,:)),0))+1])

clear cnt* corr* iClust temp unique_* Arena Control Odor PopLoc AniLoc i* curr* boot_*


% Inter-animal attraction II
% -------------------------------------------------------------------------
% (Inter-entry-Interval VS number of animals sheltered at entry)

cnt_Odor=1;
cnt_Control=1;
for iTrial = 1:length(fields(DATA.GRP.VAN.C0_01))-4
    
    % Get how many animals were sheltered at each point of time
    Odor.Pop = sum(DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Odor');
    Control.Pop = sum(DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).InShelter_Control');
    Shelter.Pop = sum([Odor.Pop; Control.Pop]);
    Arena.Pop = ones(1,DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).Duration) * DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).N - Shelter.Pop;
    
    
    OdorQuarter = [];
    ControlQuarter = [];
    for iAni = 1:DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).N
        currPos = [DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).Pos_X(:,iAni), DATA.GRP.VAN.C0_01.(['Trial_',num2str(iTrial)]).Pos_Y(:,iAni)];
        % Get location as angle
        currTheta = atan2d(currPos(:,2), currPos(:,1));
        currTheta(currTheta<0) = 360-abs(currTheta(currTheta<0));
        % Define several logicals to later sort
        OK.odor = currTheta>315 | currTheta<=45;
        OK.control = currTheta>135 & currTheta<=225;
        % Calculate the presence
        OdorQuarter = [OdorQuarter, OK.odor];
        ControlQuarter = [ControlQuarter, OK.control];
    end
    OdorQuarter = sum(OdorQuarter,2) - Odor.Pop';
    ControlQuarter = sum(ControlQuarter,2) - Control.Pop';
    clear currPos currTheta OK
    
    
    
    % Get entries
    Odor.Entries = bwlabel(diff(Odor.Pop)>0);
    Control.Entries = bwlabel(diff(Control.Pop)>0);
    
    
    % For the Odor shelter, get inter-entry-intervals and the number of
    % animals at entry
    if max(Odor.Entries) >= 2
        % Iterate over all entries
        for iEnter = 2:max(Odor.Entries)
            
            % Interval between entries
            if length(find(Odor.Entries == iEnter)) == 1
                corrOdor(2, cnt_Odor) = find(Odor.Entries == iEnter) - find(Odor.Entries == iEnter-1,1,'last');
                % Number of animals at entry
                corrOdor(1, cnt_Odor) = Odor.Pop(find(Odor.Entries == iEnter));
                % Normalize with remainnig, unsheltered animals
                corrOdor(3, cnt_Odor) = corrOdor(2, cnt_Odor) / OdorQuarter(find(Odor.Entries == iEnter)); %Arena.Pop(find(Odor.Entries == iEnter));
                % Counter
                cnt_Odor=cnt_Odor+1;
            else
                
                % Seems like animals have entered without leaving a gap.
                % Get the first interval as usual
                corrOdor(2, cnt_Odor) = find(Odor.Entries == iEnter,1) - find(Odor.Entries == iEnter-1,1,'last');
                % Number of animals at entry
                corrOdor(1, cnt_Odor) = Odor.Pop(find(Odor.Entries == iEnter,1));
                % Normalize with remainnig, unsheltered animals
                corrOdor(3, cnt_Odor) = corrOdor(2, cnt_Odor) / OdorQuarter(find(Odor.Entries == iEnter,1)); %Arena.Pop(find(Odor.Entries == iEnter,1));
                % Counter
                cnt_Odor=cnt_Odor+1;
                
                idx = find(Odor.Entries == iEnter);
                for iDouble = 2:length(idx)
                    
                    corrOdor(2, cnt_Odor) = idx(iDouble) - idx(iDouble-1);
                    % Number of animals at entry
                    corrOdor(1, cnt_Odor) = Odor.Pop(idx(iDouble));
                    % Normalize with remainnig, unsheltered animals
                    corrOdor(3, cnt_Odor) = corrOdor(2, cnt_Odor) / OdorQuarter(idx(iDouble)); %Arena.Pop(idx(iDouble));
                    
                    % Counter
                    cnt_Odor=cnt_Odor+1;
                end%iDouble
            end%if
        end%iEnter
    end%if
    
    
    % For the control shelter, get inter-entry-intervals and the number of
    % animals at entry
    if max(Control.Entries) >= 2
        % Iterate over all entries
        for iEnter = 2:max(Control.Entries)
            
            % Interval between entries
            if length(find(Control.Entries == iEnter)) == 1
                corrControl(2, cnt_Control) = find(Control.Entries == iEnter) - find(Control.Entries == iEnter-1,1,'last');
                % Number of animals at entry
                corrControl(1, cnt_Control) = Control.Pop(find(Control.Entries == iEnter));
                % Normalize with remainnig, unsheltered animals
                corrControl(3, cnt_Control) = corrControl(2, cnt_Control) / ControlQuarter(find(Control.Entries == iEnter)); %Arena.Pop(find(Control.Entries == iEnter));
                % Counter
                cnt_Control=cnt_Control+1;
            else
                
                % Seems like animals have entered without leaving a gap.
                % Get the first interval as usually
                corrControl(2, cnt_Control) = find(Control.Entries == iEnter,1) - find(Control.Entries == iEnter-1,1,'last');
                % Number of animals at entry
                corrControl(1, cnt_Control) = Control.Pop(find(Control.Entries == iEnter,1));
                % Normalize with remainnig, unsheltered animals
                corrControl(3, cnt_Control) = corrControl(2, cnt_Control) / ControlQuarter(find(Control.Entries == iEnter,1)); %Arena.Pop(find(Control.Entries == iEnter,1));
                
                % Counter
                cnt_Control=cnt_Control+1;
                
                idx = find(Control.Entries == iEnter);
                for iDouble = 2:length(idx)
                    
                    corrControl(2, cnt_Control) = idx(iDouble) - idx(iDouble-1);
                    % Number of animals at entry
                    corrControl(1, cnt_Control) = Control.Pop(idx(iDouble));
                    % Normalize with remainnig, unsheltered animals
                    corrControl(3, cnt_Control) = corrControl(2, cnt_Control) / ControlQuarter(idx(iDouble)); %Arena.Pop(idx(iDouble));
                    
                    % Counter
                    cnt_Control=cnt_Control+1;
                end%iDouble
            end%if
            
        end%iEnter
    end%if
    
end%iTrial

subplot(1,3,3)
hold on
uniquePopulation = unique(corrOdor(1,:));
cnt=1;
iAniVector = round(unique(corrOdor(1,:)),0);
for iAni = iAniVector
    
    if iAni >= 5
        idx = find(corrOdor(1,:) > iAniVector(end-3));
        QuitAfter = 1;
    else
        QuitAfter = 0;
        idx = find(corrOdor(1,:) == iAni);
    end
    
    
    currData = log10(corrOdor(2,idx)/SET.FrameRate);
    if length(currData)>1
        rng(4243)
        boot_data_odor = bootstrp(SET.N_Boot, @mean, currData);
        boot_mean = mean(boot_data_odor);
        boot_std = std(boot_data_odor);
    else
        boot_data_odor = currData;
        boot_mean = currData;
        boot_std = 0;
    end
    
    %     %--- Beeswarmplot
    %     properties.MarkerFaceColor = LocCol(find(strcmp(LocNames, 'Odor')),:);
    %     properties.MarkerSize = 2;
    %     beeswarmplot_advanced(boot_data, iAni, 0.5, properties)
    %     clear properties
    %     plot([iAni-0.25, iAni+0.25],[mean(currData), mean(currData)], 'Color', LocCol(find(strcmp(LocNames, 'Odor')),:), 'LineWidth', 2)
    
    properties.NumPoints = 250;
    properties.EdgeCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    properties.MeanWidth = 1;
    properties.OutlierCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    properties.OutlierSize = 1;
    properties.AvgType = 'mean';
    properties.MeanCol = LocCol(find(strcmp(LocNames, 'Odor')),:);
    violinplot_advanced(boot_data_odor, iAni, 0.5, properties)
    clear properties
    
    IAA_2_N.odor(cnt) = length(idx);
    cnt=cnt+1;
    
    if QuitAfter
        clear QuitAfter
        break
    end
    
end

uniquePopulation = unique(corrControl(1,:));
cnt=1;
iAniVector = round(unique(corrControl(1,:)),0);
for iAni = iAniVector
    
    if iAni >= 5
        idx = find(corrControl(1,:) > iAniVector(end-3));
        QuitAfter = 1;
    else
        QuitAfter = 0;
        idx = find(corrControl(1,:) == iAni);
    end
    
    currData = log10(corrControl(2,idx)/SET.FrameRate);
    if length(currData)>1
        rng(4243)
        boot_data_control = bootstrp(SET.N_Boot, @mean, currData);
        boot_mean = mean(boot_data_control);
        boot_std = std(boot_data_control);
    else
        boot_data_control = currData;
        boot_mean = currData;
        boot_std = 0;
    end
    
    %     %--- Beeswarmplot
    %     properties.MarkerFaceColor = LocCol(find(strcmp(LocNames, 'Control')),:);
    %     properties.MarkerSize = 2;
    %     beeswarmplot_advanced(boot_data, iAni, 0.5, properties)
    %     clear properties
    %     plot([iAni-0.25, iAni+0.25],[mean(currData), mean(currData)], 'Color', LocCol(find(strcmp(LocNames, 'Control')),:), 'LineWidth', 2)
    
    properties.NumPoints = 250;
    properties.EdgeCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    properties.MeanWidth = 1;
    properties.OutlierCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    properties.OutlierSize = 1;
    properties.AvgType = 'mean';
    properties.MeanCol = LocCol(find(strcmp(LocNames, 'Control')),:);
    violinplot_advanced(boot_data_control, iAni, 0.5, properties)
    clear properties
    
    IAA_2_N.control(cnt) = length(idx);
    cnt=cnt+1;
    
    if QuitAfter
        clear QuitAfter
        break
    end
    
end

set(gca, 'XTick', unique(corrControl(1,:)))
xlim([min(unique(corrControl(1,:)))-1, max(unique(corrControl(1,:)))+1])

clear cnt* corr* iClust temp unique_* Arena Control Odor PopLoc AniLoc i* boot_*

%% Save all figure
% Note, heatmaps were safed on the fly. we just have to create the
% colorbars

if SET.DoPlot
    % --- Colorbars
    figure
    set(gcf, 'color', 'w')
    colormap(viridis(50))
    h = colorbar;
    h.TickDirection = 'out';
    h.TickLabels = [];
    caxis([round(max(minCount),5), round(min(maxCount),5)])
    %     h.TickLabels{end} = num2str(min(maxCount));
    axis off
    export_fig('FIG\raw\cb1','-pdf')
    custom_CB = gray(1000);
    custom_CB = custom_CB(1:850,:);
    colormap(custom_CB)
    export_fig('FIG\raw\cb2','-pdf')
    
    % --- Polar_GRP_IND
    figure(hFig.Polar_GRP_IND)
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    export_fig('FIG\raw\Polar_GRP_IND','-pdf')
    
    % --- Polar_Van_FecVan
    figure(hFig.Polar_Van_FecVan)
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    export_fig('FIG\raw\Polar_Van_FecVan','-pdf')
    
    % --- Polar_all
    figure(hFig.Polar_all)
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    export_fig('FIG\raw\Polar_all','-pdf')
    
    % --- Dist2Both
    figure(hFig.Dist2Both)
    export_fig('FIG\raw\Dist2Both','-pdf')
    
    % --- PropEntry
    figure(hFig.PropEntry); hold on
    plot([0.5 8.5],[0.5 0.5],'k:')
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    xlim([0.5 8.5])
    ylim([0 1])
    export_fig('FIG\raw\PropEntry','-pdf')
    
    % --- RetentionTime
    figure(hFig.RetentionTime)
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    export_fig('FIG\raw\RetentionTime','-pdf')
    
    % --- cumTimeSheltered
    figure(hFig.cumTimeSheltered)
    subplot(1,5,1)
    ylim([0, 0.5])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [3 2 2.5 3.5])
    set(gca, 'YTick', 0:0.1:0.5)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    subplot(1,5,2)
    ylim([0, 0.5])
    xlim([0.5, 8.5])
    set(gca, 'units', 'centimeters', 'position', [9 2 1 3.5])
    set(gca, 'YTick', 0:0.1:0.5)
    set(gca, 'XTick', 1:8)
    
    subplot(1,5,3)
    ylim([-0.3, 0.3])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [15 2 2.5 3.5])
    set(gca, 'YTick', -0.3:0.1:0.3)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    subplot(1,5,4)
    ylim([-0.3, 0.3])
    xlim([0.5, 8.5])
    set(gca, 'units', 'centimeters', 'position', [21 2 1 3.5])
    set(gca, 'YTick', -0.3:0.1:0.3)
    set(gca, 'XTick', 1:8)
    
    subplot(1,5,5)
    ylim([0 1])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [27 2 3.5 3.5])
    set(gca, 'YTick', 0:0.1:1)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    export_fig('FIG\raw\cumTimeSheltered','-pdf')
    
    % --- cumTimeQuartered
    figure(hFig.cumTimeQuartered)
    subplot(1,5,1)
    ylim([0, 0.55])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [3 2 2.5 3.5])
    set(gca, 'YTick', 0:0.1:0.5)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    subplot(1,5,2)
    ylim([0, 0.55])
    xlim([0.5, 8.5])
    set(gca, 'units', 'centimeters', 'position', [9 2 1 3.5])
    set(gca, 'YTick', 0:0.1:0.5)
    set(gca, 'XTick', 1:8)
    
    subplot(1,5,3)
    ylim([-0.3, 0.3])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [15 2 2.5 3.5])
    set(gca, 'YTick', -0.3:0.1:0.3)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    subplot(1,5,4)
    ylim([-0.3, 0.3])
    xlim([0.5, 8.5])
    set(gca, 'units', 'centimeters', 'position', [21 2 1 3.5])
    set(gca, 'YTick', -0.3:0.1:0.3)
    set(gca, 'XTick', 1:8)
    
    subplot(1,5,5)
    ylim([0 1])
    xlim([0.5, SET.PropTimeBins+0.5])
    set(gca, 'units', 'centimeters', 'position', [27 2 3.5 3.5])
    set(gca, 'YTick', 0:0.1:1)
    set(gca, 'XTick', 0.5:0.5:SET.PropTimeBins+0.5)
    
    export_fig('FIG\raw\cumTimeQuartered','-pdf')
    
    % --- PropTimeSheltered
    figure(hFig.PropTimeSheltered)
    xlim([1, 17])
    ylim([0 0.5])
    set(gca, 'units', 'centimeters', 'position', [3 2 3.5 3.5])
    set(gca, 'YTick', 0:0.1:0.7, 'XTick', 2:2:16)
    
    export_fig('FIG\raw\PropTimeSheltered','-pdf')
    
    % --- Survival
    figure(hFig.Survival)
    subplot(1,3,1); hold on
    xlim([0 30])
    plot([1 30], [SET.SurvivalThreshold, SET.SurvivalThreshold], 'k:')
    set(gca, 'XScale', 'log')
    set(gca, 'YTick', 0:0.2:1)
    
    subplot(1,3,2); hold on
    xlim([0 30])
    plot([1 30], [SET.SurvivalThreshold, SET.SurvivalThreshold], 'k:')
    set(gca, 'XScale', 'log')
    set(gca, 'YTick', 0:0.2:1)
    
    subplot(1,3,3); hold on
    set(gca, 'units', 'centimeters', 'position', [50 4 3.5 3.5])
    ylim([0 20])
    export_fig('FIG\raw\Survival','-pdf')
        
    % --- Dendrogram
    figure(hFig.Dendrogram)
    set(gca, 'units', 'centimeters', 'position', [4 4 7 3.5])
    export_fig('FIG\raw\Dendrogram','-pdf')
    
    % --- SocialAttraction
    figure(hFig.SocialAttraction)
    subplot(1,3,1); hold on
    set(gca, 'units', 'centimeters', 'position', [4 4 3 4])
    subplot(1,3,2); hold on
    set(gca, 'units', 'centimeters', 'position', [10 4 3 4])
    set(gca, 'YTick', 0.5:0.25:2)
    ylim([0.5 2])
    xlim([0.5, 6.5])
    subplot(1,3,3); hold on
    set(gca, 'units', 'centimeters', 'position', [16 4 3 5])
    ylim([0 2])
    xlim([0.5, 6.5])
    set(gca, 'YTick', 0:0.5:2)
    
    export_fig('FIG\raw\SocialAttraction','-pdf')
    
    % --- PropPresence
    figure(hFig.PropPresence)
    ylim([-0.25, 0.7])
    set(gca, 'YTick', -0.2:0.1:0.75)
    xlim([1 17])
    set(gca, 'units', 'centimeters', 'position', [4 4 3.5 3.5])
    export_fig('FIG\raw\PropPresence','-pdf')
    
    % --- PreferenceIndices
    figure(hFig.PreferenceIndices)
    subplot(1,4,1)
    set(gca, 'units', 'centimeters', 'position', [4 4 3.5 3.5])
    set(gca, 'XTick', 1:8, 'YTick', -1:0.5:1)
    
    subplot(1,4,2)
    set(gca, 'units', 'centimeters', 'position', [9.5 4 3.5 3.5])
    set(gca, 'XTick', 1:8, 'YTick', -1:0.5:1)
    
    subplot(1,4,3)
    set(gca, 'units', 'centimeters', 'position', [15 4 3.5 3.5])
    set(gca, 'XTick', 1:8, 'YTick', -1:0.5:1)
    
    subplot(1,4,4)
    set(gca, 'units', 'centimeters', 'position', [20.5 4 3.5 3.5])
    set(gca, 'XTick', 1:8, 'YTick', -1:0.5:1)
    
    export_fig('FIG\raw\PreferenceIndices','-pdf')
    
    % --- Speed
    figure(hFig.Speed)
    subplot(1,2,1); hold on
    set(gca, 'units', 'centimeters', 'position', [4 4 4 7])
    subplot(1,2,2); hold on
    set(gca, 'units', 'centimeters', 'position', [10 4 4 7])
    export_fig('FIG\raw\Speed','-pdf')

end

% Save everything
save('tempDATA.mat', 'DATA', 'PooledDATA', 'AvgData', 'STATS', 'ANNOTATION', 'SET', 'maxCount', 'minCount', 'maxSurvival', 'hFig', '-v7.3')















