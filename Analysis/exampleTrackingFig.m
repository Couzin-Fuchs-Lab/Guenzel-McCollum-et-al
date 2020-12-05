%% exampleTrackingFig
% Creates examplary plots indicatiing how animals were tracked.
% 
% Verion: Final version for submission
%    30-July-2020; Yannick (MATLAB 2020a)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

%% Get data

% --- GRP tracking
GRP.Tracking = readtable('021019G12T1_tracked.csv');
GRP.Range = [17625 19125];
GRP.AniNames = unique(GRP.Tracking.id);
for iAni = 1:length(GRP.AniNames)
    idx = find(strcmp(GRP.Tracking.id, GRP.AniNames{iAni}));
    GRP.x(:,iAni) = GRP.Tracking.pos_x(idx);
    GRP.y(:,iAni) = GRP.Tracking.pos_y(idx);
end; clear iAni idx
% --- GRP frame
[handles.CurrFile, handles.CurrFilePath] = uigetfile({'*.mp4'; '*.avi'}, 'Select the GRP video file');
handles.CurrTrial = handles.CurrFile(1:end-4);
handles.Video.Obj = VideoReader([handles.CurrFilePath, handles.CurrFile]);
handles.Video.NrFrames = get(handles.Video.Obj, 'numberOfFrames')-1;
GRP.frame = uint8(rgb2gray((read(handles.Video.Obj, GRP.Range(end)))));
GRP.frame = uint8(abs(double(GRP.frame)-255));
clear handles
% --- Shelter locations
GRP.OdorLoc = [354.939161269477, 875.243499089008];
GRP.ControlLoc = [1677.79253725193, 1209.53102066197];
% --- Arena center
GRP.Center = [1025.99358084543, 1004.61778633029];
% --- Arena radius
GRP.radius = 975.3194*1.01;

% Rotate data
% --- Get the current angle between the two shelters
vec = GRP.ControlLoc - GRP.OdorLoc;
angle = atan2d(vec(2), vec(1));
M = [ cosd(-angle) -sind(-angle);
    sind(-angle)  cosd(-angle) ];
GRP.OdorLoc = (M*GRP.OdorLoc(:,1:2)')';
GRP.ControlLoc = (M*GRP.ControlLoc(:,1:2)')';
% Make sure that the odor location is right
if GRP.OdorLoc(1) < GRP.ControlLoc(1)
    disp(['Rotate GRP by ',num2str(-(180-angle)),' deg'])
else
    disp(['Rotate GRP by ',num2str(angle),' deg'])
end% if 180 rotation is needed






% --- IND tracking
IND.Tracking = readtable('040619T74_compressed_tracked.csv');
IND.Range = [1695 3195];
IND.AniNames = unique(IND.Tracking.id);
for iAni = 1:length(IND.AniNames)
    idx = find(strcmp(IND.Tracking.id, IND.AniNames{iAni}));
    IND.x(:,iAni) = IND.Tracking.pos_x(idx);
    IND.y(:,iAni) = IND.Tracking.pos_y(idx);
end; clear iAni idx
% --- IND frame
[handles.CurrFile, handles.CurrFilePath] = uigetfile({'*.mp4'; '*.avi'}, 'Select the IND video file');
handles.CurrTrial = handles.CurrFile(1:end-4);
handles.Video.Obj = VideoReader([handles.CurrFilePath, handles.CurrFile]);
handles.Video.NrFrames = get(handles.Video.Obj, 'numberOfFrames')-1;
IND.frame = uint8(rgb2gray((read(handles.Video.Obj, IND.Range(end)))));
IND.frame = uint8(abs(double(IND.frame)-255));
clear handles
% --- Shelter locations
IND.OdorLoc = [1084.19731368332, 397.789798400429];
IND.ControlLoc = [1002.62512434256, 1633.95280504451];
% --- Arena center
IND.Center = [1040.48763847060, 990.923655347685];
% --- Arena radius
IND.radius = 891.8451*1.01;

% Rotate data
% --- Get the current angle between the two shelters
vec = IND.ControlLoc - IND.OdorLoc;
angle = atan2d(vec(2), vec(1));
M = [ cosd(-angle) -sind(-angle);
    sind(-angle)  cosd(-angle) ];
IND.OdorLoc = (M*IND.OdorLoc(:,1:2)')';
IND.ControlLoc = (M*IND.ControlLoc(:,1:2)')';

% Make sure that the odor location is right
if IND.OdorLoc(1) < IND.ControlLoc(1)
    disp(['Rotate IND by ',num2str(-(180-angle)),' deg'])
else
    disp(['Rotate IND by ',num2str(angle),' deg'])
end% if 180 rotation is needed



%% Plot tracking

ColorMap = hsv(length(GRP.AniNames));

% --- GRP
figure('Name', 'GRP')
imshow(GRP.frame)
set(gcf, 'units', 'normalized', 'position', [0 0 1 1], 'color', 'w')
hold on
for iAni = 1:length(GRP.AniNames)
    plot(GRP.x(GRP.Range(1):GRP.Range(2),iAni), GRP.y(GRP.Range(1):GRP.Range(2),iAni), 'Color', ColorMap(iAni,:), 'LineWidth', 2)
    %     text(GRP.x(GRP.Range(2),iAni), GRP.y(GRP.Range(2),iAni)-14, GRP.AniNames{iAni}, 'Color', ColorMap(iAni,:), 'FontName', 'Arial', 'FontSize', 10)
end; clear iAni
axis equal
xlim([GRP.Center(1)-GRP.radius, GRP.Center(1)+GRP.radius])
ylim([GRP.Center(2)-GRP.radius, GRP.Center(2)+GRP.radius])
export_fig('GRP_tracking', '-pdf')

% --- IND
figure('Name', 'IND')
imshow(IND.frame)
set(gcf, 'units', 'normalized', 'position', [0 0 1 1], 'color', 'w')
hold on
for iAni = 1:length(IND.AniNames)
    plot(IND.x(IND.Range(1):IND.Range(2),iAni), IND.y(IND.Range(1):IND.Range(2),iAni), 'Color', ColorMap(iAni,:), 'LineWidth', 2)
    %     text(IND.x(IND.Range(2),iAni), IND.y(IND.Range(2),iAni)-14, IND.AniNames{iAni}, 'Color', ColorMap(iAni,:))
end; clear iAni
axis equal
xlim([IND.Center(1)-IND.radius, IND.Center(1)+IND.radius])
ylim([IND.Center(2)-IND.radius, IND.Center(2)+IND.radius])
export_fig('IND_tracking', '-pdf')

