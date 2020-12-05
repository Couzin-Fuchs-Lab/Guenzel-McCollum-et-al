function violinplot_advanced(ViolinData, pos, width, properties)
% VIOLINPLOT_ADVANCED(data,pos,width,properties) gives the user more
% possibilities to depict data as a violinplot. The underlying function is
% the inbuilt function ksdensity() - a kernel smoothing function estimate
% for univariate and bivariate data. It uses a probability density function
% as kernel.
%
% INPUTS:
%      - 'ViolinData':  1D vector of numbers
%      - 'pos':         x-value(vertical plot) or y-value(horizontal plot)
%                       at which the plot should be located
%      - 'width':       Total width of the violin
%      - 'properties':  Properties which specify how the violin should look 
%                       like
%
%
% DEFAULT PROPERTIES:
%       PROPERTY                        DEFAULT VALUE       DESCRIPTION
%       -----------------------------------------------------------------------------------------------------------------------------------------------------
%       properties.NumPoints =          100;                % Points at which to evaluate the probability density estimate
%       properties.MinVal =             [];                 % Smallest possible value (e.g. errors = 0)
%       properties.MaxVal =             [];                 % Biggest possible value
%       properties.Orientation =        'vertical';         %(Set how the box plot should be oriented 'vertical' or 'horizontal')
%       properties.AvgType =            'both';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
%       properties.EdgeCol =            'k';                %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
%       properties.EdgeWidth =          1;                  %(Set the width of violin's edges. If EdgeCol is 'none', this statement will be ignored)
%       properties.FillViolin =         0;                  %(Logical whether to fill the violin or not. However, this makes it impossible to save the resulting figure as a vector graphic)
%       properties.FaceCol =            'r';                %(Set the violins's face colour. Will be ignored if FillViolin is zero)
%       properties.MedianCol =          'k';                %(Set colour of median. If it should be depicted)
%       properties.MedianSymbol =       'o';                %(Set how the median should be depicted as symbol (e.g. 'x') or as line: 'line')
%       properties.MedianWidth =        5;                  %(Set width of the median's line and/or symbol)
%       properties.MeanCol =            'k';                %(Set colour of mean. If it should be depicted)
%       properties.MeanSymbol =         'line';             %(Set how the mean should be depicted as symbol (e.g. 'x') or as line: 'line')
%       properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
%       properties.SeparateOutliers =   1;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
%       properties.OutlierCol =         'k';                %(Set the colour in which outliers should be depicted)
%       properties.OutlierSymbol =      'o';                %(Set outlier symbol: e.g. 'x')
%       properties.OutlierSize =        5;                  %(Set size of outlier symbols)
%
% Version: 09-May-2020; Yannick (MATLAB 2020a)


%% Prepare data

% Replace not given property values by defaults
if nargin < 4
    properties = 'default';
end
properties = FillEmptyPropertyValues(properties);

% Exclude NaNs
data.all = ViolinData(~isnan(ViolinData));

% Get mean and median
data.median = median(data.all );
data.mean = mean(data.all );

% Get the probability density estimate based on a normal kernel
% function withequally-spaced points
[density, value] = ksdensity(data.all, 'NumPoints', properties.NumPoints);
density = density/max(density)/2*width;

% Truncate to be between min and max
density = interp1(value, density, linspace(min(data.all), max(data.all), properties.NumPoints));
value = linspace(min(data.all), max(data.all), properties.NumPoints);


% Separate outliers from data
if properties.SeparateOutliers
    
    % Define min. and max. value for data points to not be an outlier
    Q1 = quantile(data.all, 0.25);
    Q3 = quantile(data.all,0.75);
    Spread = 1.5*(Q3-Q1);
    MaxValue = Q3 + Spread;
    MinValue = Q1 - Spread;
    
    % Separate outliers from the rest of the data
    idx = find(data.all>MaxValue | data.all<MinValue);
    data.outlier = data.all(idx);
    
    % Truncate violin
    idx = find(value<MinValue | value>MaxValue);
    density(idx) = [];
    value(idx) = [];
    
    clear Q1 Q3 Spread MaxValue MinValie idx
else
    data.outlier = [];
end

% Check whether the user has set a smallest/biggest possible value
if ~isempty(properties.MinVal)
    % Truncate violin
    idx = find(value<properties.MinVal);
    density(idx) = [];
    value(idx) = [];
end
if ~isempty(properties.MaxVal)
    % Truncate violin
    idx = find(value>properties.MaxVal);
    density(idx) = [];
    value(idx) = [];
end

% Make sure we have a column vector (n,1)
density = reshape(density, [length(density),1]);
value = reshape(value, [length(value),1]);

%% Plot data
switch properties.Orientation
    case 'vertical'% Plot a vertical violinplot
       
        % Fill violin
        if properties.FillViolin
            fill([density+pos; flipud(-density+pos)], [value; flipud(value)], properties.FaceCol)
        end
        
        % Plot Violin outlines
        plot(density+pos, value, 'Color', properties.EdgeCol)
        plot(-density+pos, value, 'Color', properties.EdgeCol)
        
        
        switch properties.AvgType
            case 'median' % -----------------------------------------------                
                switch properties.MedianSymbol
                    case 'line'
                        plot([pos-width/2, pos+width/2], [data.median, data.median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
                    otherwise
                        plot(pos, data.median, properties.MedianSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
                end%MedianSymbol

            case 'mean' % -------------------------------------------------
                switch properties.MeanSymbol
                    case 'line'
                        plot([pos-width/2, pos+width/2], [data.mean, data.mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
                    otherwise
                        plot(pos, data.mean, properties.MeanSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
                end%MedianSymbol
                
            case 'both'% --------------------------------------------------
                % Median
                switch properties.MedianSymbol
                    case 'line'
                        plot([pos-width/2, pos+width/2], [data.median, data.median], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
                    otherwise
                        plot(pos, data.median, properties.MedianSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
                end%MedianSymbol
                
                % Mean
                switch properties.MeanSymbol
                    case 'line'
                        plot([pos-width/2, pos+width/2], [data.mean, data.mean], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
                    otherwise
                        plot(pos, data.mean, properties.MeanSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
                end%MedianSymbol
        end%AvgType    
        
        % Plot outlier
        plot(ones(1,length(data.outlier))*pos, data.outlier, properties.OutlierSymbol, 'MarkerFaceColor', properties.OutlierCol, 'MarkerEdgeColor', 'none', 'MarkerSize', properties.OutlierSize)

    case 'horizontal'% Plot a horizontal violinplot
                
        % Fill violin
        if properties.FillViolin
            fill([value; flipud(value)], [density+pos; flipud(-density+pos)], properties.FaceCol, 'EdgeColor', 'none')
        end
        
        % Plot Violin outlines
        plot(value, density+pos, 'Color', properties.EdgeCol)
        plot(value, -density+pos, 'Color', properties.EdgeCol)
        
        
        switch properties.AvgType
            case 'median' % -----------------------------------------------                
                switch properties.MedianSymbol
                    case 'line'
                        plot([data.median, data.median], [pos-width/2, pos+width/2], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
                    otherwise
                        plot(data.median, pos, properties.MedianSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
                end%MedianSymbol

            case 'mean' % -------------------------------------------------
                switch properties.MeanSymbol
                    case 'line'
                        plot([data.mean, data.mean], [pos-width/2, pos+width/2], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
                    otherwise
                        plot(data.mean, pos, properties.MeanSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
                end%MedianSymbol
                
            case 'both'% --------------------------------------------------
                % Median
                switch properties.MedianSymbol
                    case 'line'
                        plot([data.median, data.median], [pos-width/2, pos+width/2], 'Color', properties.MedianCol, 'LineWidth', properties.MedianWidth)
                    otherwise
                        plot(data.median, pos, properties.MedianSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MedianCol, 'MarkerSize', properties.MedianWidth)
                end%MedianSymbol
                
                % Mean
                switch properties.MeanSymbol
                    case 'line'
                        plot([data.mean, data.mean], [pos-width/2, pos+width/2], 'Color', properties.MeanCol, 'LineWidth', properties.MeanWidth)
                    otherwise
                        plot(data.mean, pos, properties.MeanSymbol, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', properties.MeanCol, 'MarkerSize', properties.MeanWidth)
                end%MedianSymbol
        end%AvgType    
        
        % Plot outlier
        plot(data.outlier, ones(1,length(data.outlier))*pos, properties.OutlierSymbol, 'MarkerFaceColor', properties.OutlierCol, 'MarkerEdgeColor', 'none', 'MarkerSize', properties.OutlierSize)
       
end%orientation



function out = FillEmptyPropertyValues(properties)
% FillEmptyPropertyValues(properties) checks whether a property has been
% set by the user. If not, replace it by its default value. If no property
% value has been set, completely use defaults.

% Set default values
default.NumPoints =          100;                % Points at which to evaluate the probability density estimate
default.MinVal =             [];                 % Smallest possible value (e.g. errors = 0)
default.MaxVal =             [];                 % Biggest possible value
default.Orientation =        'vertical';         %(Set how the box plot should be oriented 'vertical' or 'horizontal')
default.AvgType =            'both';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
default.EdgeCol =            'k';                %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
default.EdgeWidth =          1;                  %(Set the width of violin's edges. If EdgeCol is 'none', this statement will be ignored)
default.FillViolin =         0;                  %(Logical whether to fill the violin or not. However, this makes it impossible to save the resulting figure as a vector graphic)
default.FaceCol =            'r';                %(Set the violins's face colour. Will be ignored if FillViolin is zero)
default.MedianCol =          'k';                %(Set colour of median. If it should be depicted)
default.MedianSymbol =       'o';                %(Set how the median should be depicted as symbol (e.g. 'x') or as line: 'line')
default.MedianWidth =        5;                  %(Set width of the median's line and/or symbol)
default.MeanCol =            'k';                %(Set colour of mean. If it should be depicted)
default.MeanSymbol =         'line';             %(Set how the mean should be depicted as symbol (e.g. 'x') or as line: 'line')
default.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
default.SeparateOutliers =   1;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
default.OutlierCol =         'k';                %(Set the colour in which outliers should be depicted)
default.OutlierSymbol =      'o';                %(Set outlier symbol: e.g. 'x')
default.OutlierSize =        5;                  %(Set size of outlier symbols)

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
    






































