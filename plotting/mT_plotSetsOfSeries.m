function figHandle = mT_plotSetsOfSeries(PlotData, PlotStyle, varargin)
% Produces subplots plots of sets of series

% INPUT
% PlotData [subplots in y-drection by subplots in x-direction] strct array,
% with fields...
%   Xvals   Either single strcutre or [num series] long strcut array with
%           fields...
%       Vals        The i-th structure should contain the x-values for the i-th
%                   series. If the Xvals strcuture is only 1 long, these these
%                   values will be used for all series.
%   Yvals   [num series] long struct array with fields...
%       Vals
%       Fade        optional. Specifies how much to 
%       UpperError  Distance from Val
%       LowerError  Unsigned distance from Val
%   
% PlotStyle struct array with fields...
%   Xaxis [size(PlotData, 2)] struct array with fields...
%       Title       optional
%       Ticks       optional
%       TicLabels   optional
%   Yaxis [size(PlotData, 1)] struct array, with fields...
%       Title       optional
%       Ticks       optional
%       TicLabels   optional
%   Data [num series] long strcut array with fields...
%       Name        Name of the series (for legend, optional)
%       PlotType    'scatter', 'line', or 'errorShading' (shades the area
%                   in between the error bars.
%       Colour      optional  
%   Annotate [size(PlotData, 2)]*[size(PlotData, 1)] struct array with fields...
%       Text        Containing text to add to the plot
% varargin  Figure handle for the figure to plot onto. If the old figure
%           has the same subplot structure, then all the data in the old
%           subplots will be retianed.

% Joshua Calder-Travis, j.calder.travis@gmail.com


plotLineWidth = 4; %2;
axisLineWidth = 4; %2;
fontSize = 30; %18;

%% Setup

% Are we going to use the same x-values for all series?
for iSubplot = 1 : length(PlotData(:))
    
    if length(PlotData(iSubplot).Xvals) == 1
        
        PlotData(iSubplot).Xvals ...
            = repmat(PlotData(iSubplot).Xvals, length(PlotData(iSubplot).Yvals), 1);
        
    end
    
end


% Matlab uses a particualar numbering system for subplots. Find an array
% that converts from matrix index to matlab subplot number. We want an
% extra column at the end as space for a legend.
subplotHeight = size(PlotData, 1);
subplotWidth = size(PlotData, 2);

subplotIdx = NaN(subplotWidth, subplotHeight);
subplotIdx(:) = 1 : length(subplotIdx(:));
subplotIdx = subplotIdx';


% Make a new figure or use an existing one?
if isempty(varargin)
    
    figHandle = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    
else
    
    figHandle = figure(varargin{1}); 
    
end


%% Make the plots

for iPltRow = 1 : size(PlotData, 1)
    
    for iPltCol = 1 : size(PlotData, 2)
        
        % Do the plotting itself
        if (subplotHeight == 1) && (subplotWidth == 1)
            subplotObj = gca;
            hold on
            
        else
            subplotObj = subplot(subplotHeight, subplotWidth, ...
                subplotIdx(iPltRow, iPltCol));
            hold on
            
        end
        
        % Loop through all the series to be plotted
        for iSeries = 1 : length(PlotData(iPltRow, iPltCol).Yvals)
            
            if isfield(PlotStyle, 'Yaxis') ...
                    && isfield(PlotStyle.Yaxis(iPltRow), 'Ticks')

                ylim(PlotStyle.Yaxis(iPltRow).Ticks([1, end]))
                yticks(PlotStyle.Yaxis(iPltRow).Ticks)
                
            end
                
            if isfield(PlotStyle, 'Yaxis') ...
                    && isfield(PlotStyle.Yaxis(iPltRow), 'TickLabels')
                
                yticklabels(PlotStyle.Yaxis(iPltRow).TickLabels)
                
            end
            

            
            if isfield(PlotStyle, 'Xaxis') ...
                    && isfield(PlotStyle.Xaxis(iPltCol), 'Ticks')
                
                xlim(PlotStyle.Xaxis(iPltCol).Ticks([1, end]))
                xticks(PlotStyle.Xaxis(iPltCol).Ticks)
            end 
                
            if isfield(PlotStyle, 'Xaxis') ...
                    && isfield(PlotStyle.Xaxis(iPltCol), 'TickLabels')
                
                xticklabels(PlotStyle.Xaxis(iPltCol).TickLabels)
                
            end
            
            if isfield(PlotStyle, 'Xaxis') ...
                    && isfield(PlotStyle.Xaxis(iPltCol), 'Title')
                
                xLabel = PlotStyle.Xaxis(iPltCol).Title;
                
            else
                xLabel = [];
                
            end

            
            % Whether we want to add labels depends on whether we are at the
            % edge of the figure
            if isfield(PlotStyle, 'Yaxis') ...
                    && isfield(PlotStyle.Yaxis(iPltRow), 'Title')
                
                yLabel = PlotStyle.Yaxis(iPltRow).Title;
                
            else
                yLabel = [];
                
            end
            
            
            if iPltCol == 1; ylabel(yLabel, 'FontSize',fontSize); end
            if iPltRow == size(PlotData, 1); xlabel(xLabel, 'FontSize',fontSize); end
            

            % What colour should we plot in?
            if isfield(PlotStyle.Data, 'Colour') ...
                && ~isempty(PlotStyle.Data(iSeries).Colour)
                
                plottingColour = PlotStyle.Data(iSeries).Colour;
                
            % Otherwise use default colours
            else
                for iColour = 1 : 5
                    defaultColours{iColour}  = mT_pickColour(iColour);
                end
                
                plottingColour = defaultColours{iSeries};
                
                % Store for legend
                PlotStyle.Data(iSeries).Colour = plottingColour;
                
            end

            
            
            
            % Has fading been requested? If not set all fade values to 1
            if ~isfield(PlotStyle, 'Yaxis') ...
                    || ~isfield(PlotData(iPltRow, iPltCol).Yvals(iSeries), ...
                    'Fade') ...
                    || isempty(PlotData(iPltRow, iPltCol).Yvals(iSeries).Fade)
                
                PlotData(iPltRow, iPltCol).Yvals(iSeries).Fade = ...
                    ones(size( ...
                    PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals));
                
                
            end
            
            
            % Loop through and plot every point
            for iXpos = 1 : length(PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals)
                
                % What is the fading for this point. Note, not used for
                % plot type, scatter.
                colourIncFade = plottingColour;
                colourIncFade(end +1) = ...
                    PlotData(iPltRow, iPltCol).Yvals(iSeries).Fade(iXpos);
                
                currentX = ...
                    PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals(iXpos);
                currentY = ...
                    PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos);
                    
                
                % What type of plot are we doing?
                if strcmp(PlotStyle.Data(iSeries).PlotType, 'scatter')
                    
                    scatter(currentX, currentY, 'o', ...
                        'MarkerEdgeColor', plottingColour, ...
                        'LineWidth', plotLineWidth)
                    
                    errorbar(currentX, currentY, ...
                        PlotData(iPltRow, iPltCol...
                            ).Yvals(iSeries).UpperError(iXpos), ...
                        PlotData(iPltRow, iPltCol...
                            ).Yvals(iSeries).LowerError(iXpos), ...
                        'LineStyle','none', 'Color', plottingColour, ...
                        'LineWidth', plotLineWidth);
                    
                
                elseif strcmp(PlotStyle.Data(iSeries).PlotType, 'line')
                % Fading will determine half of the line/shading to
                % the right and half of the line/shading to the left of the
                % Xval, for line plots.    
                
                currentX = ...
                    PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals(iXpos);
                currentY = ...
                    PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos);
                        
                
                    % Start with the lines to the right of the point. There
                    % is nothing to draw to the right if this is the final
                    % data point.
                    if iXpos < ...
                            length(PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals)
                        
                        nextX = ...
                            PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals(iXpos +1);
                        nextY = ...
                            PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos +1);
                        
                        plot([currentX, (currentX + nextX)/2], ...
                            [currentY, (currentY + nextY)/2], ...
                            'Color', colourIncFade, ...
                            'LineWidth', plotLineWidth);
                    
                    end


                    % Now the lines to the left of the point. There
                    % is nothing to draw to the right if this is the final
                    % data point.
                    if iXpos > 1
                        
                        prevX = ...
                            PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals(iXpos -1);
                        prevY = ...
                            PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos -1); 
                        
                        plot([(prevX + currentX)/2, currentX], ...
                            [(prevY + currentY)/2, currentY], ...
                            'Color', colourIncFade, ...
                            'LineWidth', plotLineWidth);
                    
                    end
                    
                    
                elseif strcmp(PlotStyle.Data(iSeries).PlotType, ...
                        'errorShading')
                    
                    currentError = [PlotData(iPltRow, iPltCol ...
                            ).Yvals(iSeries).LowerError(iXpos), ...
                        -PlotData(iPltRow, iPltCol ...
                            ).Yvals(iSeries).UpperError(iXpos)] + ...
                            currentY;
                         
                        
                    if iXpos < length(PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals)
                        
                        nextX = ...
                            PlotData(iPltRow, iPltCol).Xvals(iSeries).Vals(iXpos +1);
                        nextY = ...
                            PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos +1);
                        
                        
                        nextError = [PlotData(iPltRow, iPltCol ...
                            ).Yvals(iSeries).LowerError(iXpos +1), ...
                        -PlotData(iPltRow, iPltCol ...
                            ).Yvals(iSeries).UpperError(iXpos +1)] + ...
                            nextY;

                        fill([currentX, currentX, nextX, nextX], ...
                            [currentError, fliplr(nextError)], ...
                            plottingColour, 'Edgecolor','none', ...
                            'FaceAlpha',.25);
                    
                    end


%                     % Now the shading to the left of the point. There
%                     % is nothing to draw to the right if this is the final
%                     % data point.
%                     if iXpos > 1
%                         
%                         prevX = PlotData(iPltRow, iPltCol).Xvals(iXpos -1);
%                         prevY = ...
%                             PlotData(iPltRow, iPltCol).Yvals(iSeries).Vals(iXpos -1); 
%                         
%                         
%                         prevError = [PlotData(iPltRow, iPltCol ...
%                             ).Yvals(iSeries).LowerError(iXpos -1), ...
%                         -PlotData(iPltRow, iPltCol ...
%                             ).Yvals(iSeries).UpperError(iXpos -1)] + ...
%                             prevY;
%           
%                         fill([prevX, prevX, nextX, nextX], ...
%                             [prevError, fliplr(currentError)], colourIncFade);
%                     
%                     end
                     
                end
                
            end
            
            subplotObj.LineWidth = axisLineWidth;
            subplotObj.FontSize = fontSize;
            
        end
        
        
        % Add annotations
        if isfield(PlotStyle, 'Annotate')
        
        loc = get(subplotObj, 'position');
        dim = loc .* [1 1 0.25 0.25];
%         xlimits = xlim;
%         ylimits = ylim;
% [xlimits(1), ylimits(2) - ((ylimits(2)-ylimits(1))/4), ...
%             (xlimits(2)-xlimits(1))/4, ((ylimits(2)-ylimits(1))/4)]
%         
        annotation('textbox', dim, ...
            'String', PlotStyle.Annotate(iPltRow, iPltCol).Text, ...
            'LineStyle', 'none')
        
        end
        
    end
    
end
                    

%% Make the legends

% Making the legend is a little complicated in this case. We are going
% to trick MATLAB and draw invisible new lines.

% Set legend
if isfield(PlotStyle, 'Data') && isfield(PlotStyle.Data, 'Name')

    if (subplotHeight == 1) && (subplotWidth == 1)
        % Do nothing, there is only one plot to pick from
        
    else
        subplot(subplotHeight, subplotWidth, subplotIdx(ceil(subplotHeight/2), end));
        hold on
    end
    
    
    legendLabels = cell(1, length(PlotStyle.Data));
    
    for iLabel = 1 : length(legendLabels)
        
        legendLabels{iLabel} = PlotStyle.Data(iLabel).Name;
        
        legendLine(iLabel) = ...
            plot(NaN, NaN, 'Marker', 'o', 'Color', PlotStyle.Data(iLabel).Colour, ...
            'LineWidth', plotLineWidth*3);
        
        
    end
    
    
    legObj = legend(legendLine, legendLabels{:});
    legend boxoff
    
    legObj.FontSize = fontSize;
    legObj.LineWidth = axisLineWidth;


end

% title(legObj, 'Display items')
% legObj.Title.Visible = 'on';
% legObj.FontSize = 14;


