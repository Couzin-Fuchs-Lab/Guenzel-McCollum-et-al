function boxplot_advanced(data, pos, width, properties)
% BOXPLOT_ADVANCED(data,pos,width,properties) gives the user more
% possibilities to depict data as a boxplot. For a detailed description of 
% possible propperties, together with example values, see the list below.
% - Whiskers show minimum and maximum data points (outliers excluded).
% - Notches represent the 95% confidence interval of the median. The notch
%   width (i.e. how much they drag into the box) is set the 1/4 of the
%   total box width per notch
%
% INPUTS:
%      - 'data':        1D vector of numbers
%      - 'pos':         x-value(vertical plot) or y-value(horizontal plot)
%                       at which the plot should be located
%      - 'width':       Total width of the box
%      - 'properties':  Properties which specify how the boxplot should look like
%
% 
% DEFAULT PROPERTIES:
%       PROPERTY                        DEFAULT VALUE       DESCRIPTION
%       -----------------------------------------------------------------------------------------------------------------------------------------------------
%       properties.Orientation =        'vertical';         %(Set how the box plot should be oriented'vertical' or 'horizontal')
%       properties.Type =               'median';           %(Set which measure should be plotted: 'median', 'mean' or 'both')
%       properties.Notch =              'false';            %(Set whether to plot notches: 'false', or 'true')
%       properties.NotchWidth =         0.25;               %(Set notch width (per notch) as a fraction of the total box width. Maximum possible value is 0.5)
%       properties.BoxEdgeCol =         'none';             %(Color classification of the box's edge. Set to 'none' if no edge should be depicted)
%       properties.BoxEdgeWidth =       1;                  %(Set the width of box's edges. If BoxEdgeCol is 'none', this statement will be ignored)
%       properties.BoxFaceCol =         'k';                %(Set the box's face colour)
%       properties.MedianCol =          'w';                %(Set colour of median. If it should be depicted)
%       properties.MedianSymbol =       'line';             %(Set how the median should be depicted: as a line: 'line', as a symbol: 'x', or both 'line+x')
%       properties.MedianWidth =        2;                  %(Set width of the median's line and/or symbol)
%       properties.MeanCol =            'w';                %(Set colour of mean. If it should be depicted)
%       properties.MeanSymbol =         'line';             %(Set how the mean should be depicted: as a line: 'line', as a symbol: 'x', or both 'line+x')
%       properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
%       properties.OutlierCol =         'k';                %(Set the colour in which outliers should be depicted)
%       properties.OutlierSymbol =      '+';                %(Set outlier symbol: 'x')
%       properties.OutlierSize =        12;                 %(Set size of outlier symbols)
%       properties.OutlierWidth =       1.5;                %(Set width of outlier symbols)
%       properties.WhiskerStyle =       '-';                %(Set style of line for max and min from box)
%       properties.WhiskerWidth =       2;                  %(Set the sdjacents' widths)
%       properties.WhiskerCol =         'k';                %(Set color of Whiskers)
%       properties.WhiskerTopLength =   0.6*width;          %(Set length of the Whisker top)
%
% -------------------------------------------------------------------------
% HOW TO READ RESULTING PLOTS
% 
%           +                   Outlier = value bigger than 1.5*IQR
%
%
%      ___________              Maximum value (outliers excluded)
%           |
%           |
%    _______|_______            Q3
%   |               |
%   |               |
%   |               |           95% CI of the median (upper margin)
%    \             /            
%     \           /             
%      \         /
%       ----o----               Median (depicted as line and symbol
%      /         \              
%     /           \             
%    /             \
%   |               |           95% CI of the median (lower margin)
%   |               |
%   |               |
%   |_______________|           Q1
%           |
%           |
%           |
%      _____|_____              Minimum value (outliers excluded)
%
%
%           +                   Outlier = value smaller than 1.5*IQR
%
% Version: Yannick, 18-June-18

% Replace not given property values by defaults
if nargin < 4
    properties = 'default';
end
properties = FillEmptyPropertyValues(properties, width);

% Tell the user, that if notches are enabled, the line indicating the 
% mean's position might exeed the box's edges. 
if (strcmp(properties.Type, 'mean') || strcmp(properties.Type, 'both')) && length(properties.MedianSymbol)<=4
        warning('Line indicating position of the data"s mean might exceed the box?s edges!')
end

% Exclude NaNs
data(isnan(data)) = [];

% Calculate values of the boxplot (Median, quantiles, spread, MinValue and MaxValue)
Median = median(data);
Mean = mean(data);
Q1 = quantile(data, 0.25);
Q3 = quantile(data,0.75);
Spread = 1.5*(Q3-Q1);
MaxValue = Q3 + Spread;
MinValue = Q1 - Spread;

% Get outliers
Outliers = data( data>MaxValue | data<MinValue);

% Set length of the median
median_length = [pos-(width/2), pos+(width/2)];
if strcmp(properties.Notch, 'true')
    if properties.NotchWidth > 0.5
        properties.NotchWidth = 0.5;
        warning('The notch width exeeds half of the box?s width! Thus, it is set to its maximum possible value of 0.5')
    end%if notch width too big
    median_length = [pos-(width/2-width*properties.NotchWidth), pos+(width/2-width*properties.NotchWidth)];    
end%if

% Plot
hold on
%--------------------------------------------------------------------------
if strcmp(properties.Orientation, 'vertical')% Plot a vertical boxplot
    
    % If there are outliers, plot them
    if~isempty(Outliers)
        x = ones(1,length(Outliers))*pos;
        for i = 1:length(x)
            plot(x(i), Outliers(i), 'Color', properties.OutlierCol, 'Marker', properties.OutlierSymbol, 'MarkerSize', properties.OutlierSize, 'LineWidth', properties.OutlierWidth)
        end        
        % Adjust max and min values
        MaxValue = max(setdiff(data, Outliers));
        MinValue = min(setdiff(data, Outliers));
    elseif isempty(Outliers)
        % Adjust max and min values
        MaxValue = max(data);
        MinValue = min(data);
    end
    
    % Plot upper and lower Whisker
    plot([pos, pos], [Q3, MaxValue], 'Color', properties.WhiskerCol, 'LineStyle', properties.WhiskerStyle, 'LineWidth',properties.WhiskerWidth)
    plot([pos, pos], [Q1, MinValue], 'Color', properties.WhiskerCol, 'LineStyle', properties.WhiskerStyle, 'LineWidth',properties.WhiskerWidth)
    plot([pos-(properties.WhiskerTopLength/2), pos+(properties.WhiskerTopLength/2)], [MaxValue, MaxValue], 'Color', properties.WhiskerCol, 'LineWidth',properties.WhiskerWidth)
    plot([pos-(properties.WhiskerTopLength/2), pos+(properties.WhiskerTopLength/2)], [MinValue, MinValue], 'Color', properties.WhiskerCol, 'LineWidth',properties.WhiskerWidth)
    
    % Plot box
    PlotBox(properties, pos, width, data, Q1, Q3, Median, median_length)
    
    % Add an option to choose between mean, median or both to plot 
    if strcmp(properties.Type, 'median')
        
        % - - - Plot Median - - -
        if strcmp(properties.MedianSymbol, 'line')
            % As a line
            plot(median_length, [Median, Median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)            
        elseif length(properties.MedianSymbol) >5 && strcmp(properties.MedianSymbol(1:5), 'line+')
            % As line and symbol
            plot(median_length, [Median, Median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
            plot(pos, Median, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol(end), 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)            
        else
            % Only as a symbol
            plot(pos, Median, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol, 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
        end
        
    elseif strcmp(properties.Type, 'mean') %-------------------------------
        
        % - - - Plot Mean - - -
        if strcmp(properties.MeanSymbol, 'line')
            % As a line
            plot([pos-(width/2), pos+(width/2)],[Mean, Mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)            
        elseif length(properties.MeanSymbol) >5 && strcmp(properties.MeanSymbol(1:5), 'line+')
            % As line and symbol
            plot([pos-(width/2), pos+(width/2)],[Mean, Mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
            plot(pos, Mean, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol(end), 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)            
        else
            % Only as a symbol
            plot(pos, Mean, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol, 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
        end
        
    elseif strcmp(properties.Type, 'both') %-------------------------------
        
        % First, plot the median (same code as above for median only, i.e. same properties
        if strcmp(properties.MedianSymbol, 'line')
            % As a line
            plot(median_length, [Median, Median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)            
        elseif length(properties.MedianSymbol) >5 && strcmp(properties.MedianSymbol(1:5), 'line+')
            % As line and symbol
            plot(median_length, [Median, Median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
            plot(pos, Median, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol(end), 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)            
        else
            % Only as a symbol
            plot(pos, Median, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol, 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
        end
        
        % Second, plot the mean (same code as above for mean only, i.e. same properties
        if strcmp(properties.MeanSymbol, 'line')
            % As a line
            plot([pos-(width/2), pos+(width/2)],[Mean, Mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)            
        elseif length(properties.MeanSymbol) >5 && strcmp(properties.MeanSymbol(1:5), 'line+')
            % As line and symbol
            plot([pos-(width/2), pos+(width/2)],[Mean, Mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
            plot(pos, Mean, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol(end), 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)            
        else
            % Only as a symbol
            plot(pos, Mean, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol, 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
        end%if     
        
    end% if properties.Type
    
%--------------------------------------------------------------------------
elseif strcmp(properties.Orientation, 'horizontal')% Plot a horizontal boxplot
    
    % If there are outliers, plot them
    if~isempty(Outliers)
        y = ones(1,length(Outliers))*pos;
        for i = 1:length(y)
            plot(Outliers(i), y(i), 'Color', properties.OutlierCol, 'Marker', properties.OutlierSymbol, 'MarkerSize', properties.OutlierSize, 'LineWidth', properties.OutlierWidth)
        end
        % Adjust max and min values
        MaxValue = max(setdiff(data, Outliers));
        MinValue = min(setdiff(data, Outliers));
    elseif isempty(Outliers)
        % Adjust max and min values
        MaxValue = max(data);
        MinValue = min(data);
    end%if
    
    % Plot upper and lower Whisker
    plot([Q3, MaxValue], [pos, pos], 'Color', properties.WhiskerCol, 'LineStyle', properties.WhiskerStyle, 'LineWidth',properties.WhiskerWidth)
    plot([Q1, MinValue], [pos, pos], 'Color', properties.WhiskerCol, 'LineStyle', properties.WhiskerStyle, 'LineWidth',properties.WhiskerWidth)
    plot([MaxValue, MaxValue], [pos-(properties.WhiskerTopLength/2), pos+(properties.WhiskerTopLength/2)], 'Color', properties.WhiskerCol, 'LineWidth',properties.WhiskerWidth)
    plot([MinValue, MinValue], [pos-(properties.WhiskerTopLength/2), pos+(properties.WhiskerTopLength/2)], 'Color', properties.WhiskerCol, 'LineWidth',properties.WhiskerWidth)
    
    % Plot box
    PlotBox(properties, pos, width, data, Q1, Q3, Median, median_length)
    
    if strcmp(properties.Type, 'median')
        
        % - - - Plot Median - - -
        if strcmp(properties.MedianSymbol, 'line')
            % As a line
            plot([Median, Median], median_length, 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)            
        elseif length(properties.MedianSymbol) >5 && strcmp(properties.MedianSymbol(1:5), 'line+')
            % As line and symbol
            plot([Median, Median], median_length, 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
            plot(Median, pos, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol(end), 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
        else
            % Only as a symbol
            plot(Median, pos, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol, 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
        end
        
    elseif strcmp(properties.Type, 'mean') %-------------------------------
        
        % - - - Plot Mean - - -
        if strcmp(properties.MeanSymbol, 'line')
            % As a line
            plot([Mean, Mean], [pos-(width/2), pos+(width/2)], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)            
        elseif length(properties.MeanSymbol) >5 && strcmp(properties.MeanSymbol(1:5), 'line+')
            % As line and symbol
            plot([Mean, Mean], [pos-(width/2), pos+(width/2)], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
            plot(Mean, pos, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol(end), 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)    
        else
            % Only as a symbol
            plot(Mean, pos, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol, 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
        end
        
    elseif strcmp(properties.Type, 'both') %-------------------------------
        
        % First, plot the median (same code as above for median only, i.e. same properties
        if strcmp(properties.MedianSymbol, 'line')
            % As a line
            plot([Median, Median], median_length, 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)            
        elseif length(properties.MedianSymbol) >5 && strcmp(properties.MedianSymbol(1:5), 'line+')
            % As line and symbol
            plot([Median, Median], median_length, 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
            plot(Median, pos, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol(end), 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)  
        else
            % Only as a symbol
            plot(Median, pos, 'Color', properties.MedianCol, 'Marker', properties.MedianSymbol, 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
        end%if
        
        % Second, plot the mean (same code as above for mean only, i.e. same properties
        if strcmp(properties.MeanSymbol, 'line')
            % As a line
            plot([Mean, Mean], [pos-(width/2), pos+(width/2)], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)            
        elseif length(properties.MeanSymbol) >5 && strcmp(properties.MeanSymbol(1:5), 'line+')
            % As line and symbol
            plot([Mean, Mean], [pos-(width/2), pos+(width/2)], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
            plot(Mean, pos, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol(end), 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
        else
            % Only as a symbol
            plot(Mean, pos, 'Color', properties.MeanCol, 'Marker', properties.MeanSymbol, 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
        end%if
        
    end%if properties.Type
    
end%if vertical or horizontal
clear

function out = FillEmptyPropertyValues(properties, width)
% FillEmptyPropertyValues(properties) checks whether a property has been
% set by the user. If not, replace it by its default value. If no property
% value has been set, completely use defaults.

% Set default values
default.Orientation =        'vertical';         %(Set how the box plot should be oriented'vertical' or 'horizontal')
default.Type =               'median';           %(Set which measure should be plotted: 'median', 'mean' or 'both')
default.Notch =              'false';            %(Set whether to plot notches: 'false', or 'true')
default.NotchWidth =         0.25;               %(Set notch width (per notch) as a fraction of the total box width. Maximum possible value is 0.5)
default.BoxEdgeCol =         'none';             %(Color classification of the box's edge. Set to 'none' if no edge should be depicted)
default.BoxEdgeWidth =       1;                  %(Set the width of box's edges. If BoxEdgeCol is 'none', this statement will be ignored)
default.BoxFaceCol =         'k';                %(Set the box's face colour)
default.MedianCol =          'w';                %(Set colour of median. If it should be depicted)
default.MedianSymbol =       'line';             %(Set how the median should be depicted: as a line: 'line', as a symbol: 'x', or both 'line+x')
default.MedianWidth =        2;                  %(Set width of the median's line and/or symbol)
default.MeanCol =            'w';                %(Set colour of mean. If it should be depicted)
default.MeanSymbol =         'line';             %(Set how the mean should be depicted: as a line: 'line', as a symbol: 'x', or both 'line+x')
default.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
default.OutlierCol =         'k';                %(Set the colour in which outliers should be depicted)
default.OutlierSymbol =      '+';                %(Set outlier symbol: 'x')
default.OutlierSize =        12;                 %(Set size of outlier symbols)
default.OutlierWidth =       1.5;                %(Set width of outlier symbols)
default.WhiskerStyle =       '-';                %(Set style of line for max and min from box)
default.WhiskerWidth =       2;                  %(Set the sdjacents' widths)
default.WhiskerCol =         'k';                %(Set color of Whiskers)
default.WhiskerTopLength =   0.6*width;         %(Set length of the Whisker top)

% If the user has not given any property, completely use the default
% settings. Alternatively, check for missing property values and fill gaps
% with default values
if strcmp(properties, 'default')
    out = default;
else
    % Get field names in order to be able to iterate over each property to
    % check whether it has been given by the user
    FieldNames = fieldnames(default);
    
    % Iterate over all field names and check whether to fill a gap or not
    for iProp = 1:length(FieldNames)
        if ~isfield(properties, FieldNames{iProp})
            out.(FieldNames{iProp}) = default.(FieldNames{iProp});
        else
            out.(FieldNames{iProp}) = properties.(FieldNames{iProp});
        end%if property has not been set
    end%iProp
    
end%if completely use default values


function PlotBox(properties, pos, width, data, Q1, Q3, Median, median_length)
% PlotBox(Median, Q1, Q3, properties) plots the box with or without notches


% Either plot a vertical or horizontal box with or without notches
if strcmp(properties.Notch, 'false') %-------------------------------------
    
    % Vertical, w/o notches -----------------------------------------------
    if strcmp(properties.Orientation, 'vertical')
        x = pos-(width/2);
        y = Q1;
        h = Q3-Q1;
        rectangle('Position', [x,y,width,h], 'FaceColor', properties.BoxFaceCol, 'EdgeColor', properties.BoxEdgeCol, 'LineWidth',properties.BoxEdgeWidth)
    
    % Horizontal, w/o notches ---------------------------------------------
    elseif strcmp(properties.Orientation, 'horizontal')
        x = Q1;
        y = pos-(width/2);
        w = Q3-Q1;
        h = width;
        rectangle('Position', [x,y,w,h], 'FaceColor', properties.BoxFaceCol, 'EdgeColor', properties.BoxEdgeCol, 'LineWidth',properties.BoxEdgeWidth)
    end
    
elseif strcmp(properties.Notch, 'true')% ----------------------------------
    
    CI = ConfidenceInterval_median(data);

    % Vertical, w/ notches ------------------------------------------------
    if strcmp(properties.Orientation, 'vertical')
        X = [pos-(width/2), pos+(width/2), pos+(width/2), median_length(2), pos+(width/2), pos+(width/2), pos-(width/2), pos-(width/2), median_length(1), pos-(width/2)];
        Y = [Q3, Q3, CI(2), Median, CI(1), Q1, Q1, CI(1), Median, CI(2)];
        h = patch(X, Y, 'k');
        h.FaceColor = properties.BoxFaceCol;
        h.EdgeColor = properties.BoxEdgeCol;
        h.LineWidth = properties.BoxEdgeWidth;
    
    % Horizontal, w/ notches ----------------------------------------------
    elseif strcmp(properties.Orientation, 'horizontal')
        X = [Q1, CI(1), Median, CI(2), Q3, Q3, CI(2), Median, CI(1), Q1];
        Y = [pos+(width/2), pos+(width/2), median_length(2), pos+(width/2), pos+(width/2), pos-(width/2), pos-(width/2), median_length(1), pos-(width/2), pos-(width/2)];
        h = patch(X, Y, 'k');
        h.FaceColor = properties.BoxFaceCol;
        h.EdgeColor = properties.BoxEdgeCol;
        h.LineWidth = properties.BoxEdgeWidth;
    end
    
end


function CI = ConfidenceInterval_median(data)
% Compute the 95% confidence interval (CI) on the median of a 1D vector

% Set Z for a 95% CI
Z = 1.96;
% Preallocation of later output
CI = NaN(2,1);
% Sort data
temp_data = sort(data);
% Calculate upper and lower CI margin
CI(1) = temp_data(round((length(data)/2) - ((Z*sqrt(length(data)))/2)));
CI(2) = temp_data(round(1+ (length(data)/2) + ((Z*sqrt(length(data)))/2)));