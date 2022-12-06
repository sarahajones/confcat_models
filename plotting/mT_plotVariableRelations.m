function figHandle = mT_plotVariableRelations(DSet, XVars, YVars, Series, ...
    PlotStyle, varargin)
% Plots various every relationship between the x-variables variables 
% and the y-variables (after binning the x-variables) in subplots.


% INPUT
% XVars     [num x-variables] long struct array with fields...
%   ProduceVar      Function handle. Function accepts 'DSet.P(i).Data' (see
%                   standard data strcuture information in repo README).
%                   Function produces [num trials] long vector of values of this
%                   variable.
%   NumBins         When the x-variable is binned for plotting and calculation
%                   of the y-variable, how many bins should be used? If the data
%                   is already binned, pass 'prebinned'.
%   FindIncludedTrials
%                   (optional)
%                   Function handle. Function accepts 'DSet.P(i).Data' and 
%                   produces a [num trials] long logical of trials to include
%                   when evaluating this x-variable. Note that trials are only
%                   included if they meet all inclusion criteria, see 'Series'.
% 
% YVars     [num y-variables] long struct array with fields...
%   ProduceVar      Function handle. Function accepts 'DSet.P(i).Data' (see
%                   standard data strcuture information in repo README), and 
%                   'binTrials' a [num trials] long logical, with a one for
%                   every trial in the bin currently being evaluated, that 
%                   meets all inclusions criteria (for the YVar and the Series).
%                   Function produces a single value, the value of the
%                   y-variable for this bin.
%   FindIncludedTrials
%                   Function handle. Function accepts 'DSet.P(i).Data' and 
%                   produces a [num trials] long logical of trials to include
%                   when evaluating this y-variable. Note that trials are only
%                   included if they meet all inclusion criteria, see 'Series'.
%
% Series    [num series] long struct array. Each series will be plotted in 
%           every subplot. Contains fields...
%   FindIncludedTrials
%                   Function handle. Function accepts 'DSet.P(i).Data' and 
%                   produces a [num trials] long logical of trials to include
%                   when evaluating this series. Note that trials are only
%                   included if they meet all inclusion criteria, see 'YVars'.
%
% PlotStyle struct array with fields...
%   Xaxis [size(PlotData, 2)] struct array with fields...
%       Title       optional
%       Ticks       optional
%       TickLabels   optional
%   Yaxis %   Yaxis [size(PlotData, 1)] struct array with fields...
%       Title       optional
%       Ticks       optional
%       TickLabels   optional
%   Data [num series] long strcut array with fields...
%       Name        Name of the series (for legend, optional)
%       PlotType    'scatter', 'line', or 'errorShading' (shades the area
%                   in between the error bars.
%       Colour        
% varargin  Figure handle for the figure to plot onto. If the old figure
%           has the same subplot structure, then all the data in the old
%           subplots will be retianed.

% Joshua Calder-Travis, j.calder.travis@gmail.com


if ~isempty(varargin)
    figHandle = varargin{1};
    hold on
else
    figHandle = figure;
end


%% Binning per participant

% Initialise
PtpntPlotData.BinX = [];
PtpntPlotData.BinY = [];
PtpntPlotData = repmat(PtpntPlotData, length(YVars), length(XVars));

for iY = 1 : length(YVars)
    
    for iX = 1 : length(XVars)
        
        % How many bins to we have/want?
        if ~strcmp(XVars(iX).NumBins, 'prebinned')
            
            numBins = XVars(iX).NumBins;
            
        elseif strcmp(XVars(iX).NumBins, 'prebinned')
            
            numBins = determineBinNumber(DSet, XVars, YVars, Series, iX);
            
        end

        PtpntPlotData(iY, iX).BinX = ...
            NaN(length(Series), numBins, length(DSet.P));
        PtpntPlotData(iY, iX).BinY = ...
            NaN(length(Series), numBins, length(DSet.P));
        
    end
    
end


for iPtpnt = 1 : length(DSet.P)
    
    for iY = 1 : length(YVars)
        
        for iS = 1 : length(Series)
            
            for iX= 1 : length(XVars)
                
                xValues = XVars(iX).ProduceVar(DSet.P(iPtpnt).Data);
                
                % Which data meets all inclusion conditions?
                includedData = findIncludedTrials(DSet, iPtpnt, XVars, iX, YVars, iY, Series, iS);
                
                % Bin data meeting the inclusion conditions
                if ~strcmp(XVars(iX).NumBins, 'prebinned')
                    
                    % To exclude data we pass the 'blockType' argument of the below
                    % function as NaNs.
                    blockType = double(includedData);
                    blockType(blockType == 0) = NaN;
                    
                    % Is there any data to bin?
                    if sum(~isnan(blockType))==0
                        warning(['No data for x ' num2str(iX) ...
                            ', y ' num2str(iY) ' s, ' num2str(iS) '.'])
                        continue
                    end
                    
                    
                    BinSettings.DataType = 'integer';
                    BinSettings.BreakTies = false;
                    BinSettings.Flip = false;
                    BinSettings.EnforceZeroPoint = false;
                    BinSettings.NumBins = XVars(iX).NumBins;
                    BinSettings.SepBinning = false;
                    
                    [ordinalVar, ~, ~] = mT_makeVarOrdinal(BinSettings, xValues, ...
                        blockType, []);
                    
                    
                else
                    
                    ordinalVar = xValues;
                    ordinalVar(~includedData) = NaN;
                    
                end
                
                
                bins = unique(ordinalVar);
                bins(isnan(bins)) = [];
                
                for iBin = 1 : length(bins)
                    
                    binTrials = ordinalVar == bins(iBin);
                    
                    % We will use the average xVar value of the data points in
                    % the bin as the bin's x-position
                    PtpntPlotData(iY, iX).BinX(iS, iBin, iPtpnt) ...
                        = mean(xValues(binTrials));
                    
                    % Compute the y-position
                    PtpntPlotData(iY, iX).BinY(iS, iBin, iPtpnt) ...
                        = YVars(iY).ProduceVar(DSet.P(iPtpnt).Data, binTrials);
                    
                end
                
            end
            
        end
        
    end
    
end
      

%% Averaging over participants

% We now average the results over participants, taking the mean bin x- and
% y-positions.
for iY = 1 : length(YVars)
    
    for iX = 1 : length(XVars)
        
        averageXdata = mean(PtpntPlotData(iY, iX).BinX(:, :, :), 3);
        averageYdata = mean(PtpntPlotData(iY, iX).BinY(:, :, :), 3);
        SEM = std(PtpntPlotData(iY, iX).BinY(:, :, :), 0, 3) ...
            ./ (sum(~isnan(PtpntPlotData(iY, iX).BinY(:, :, :)), 3).^(1/2));
        
        for iS = 1 : length(Series)
        
            AvPlotData(iY, iX).Xvals(iS).Vals = averageXdata(iS, :);
            AvPlotData(iY, iX).Yvals(iS).Vals = averageYdata(iS, :);
            AvPlotData(iY, iX).Yvals(iS).UpperError = SEM(iS, :);
            AvPlotData(iY, iX).Yvals(iS).LowerError = SEM(iS, :);
            
        end
        
    end
    
end


figHandle = mT_plotSetsOfSeries(AvPlotData, PlotStyle, figHandle);  


end


function includedData = findIncludedTrials(DSet, iPtpnt, XVars, iX, YVars, iY, Series, iS)
% Find the trials meeting all inclusion criteria

if isfield(XVars, 'FindIncludedTrials')
    
    includedData ...
        = XVars(iX).FindIncludedTrials(DSet.P(iPtpnt).Data) ...
        & YVars(iY).FindIncludedTrials(DSet.P(iPtpnt).Data) ...
        & Series(iS).FindIncludedTrials(DSet.P(iPtpnt).Data);
    
else
    
    includedData ...
        = YVars(iY).FindIncludedTrials(DSet.P(iPtpnt).Data) ...
        & Series(iS).FindIncludedTrials(DSet.P(iPtpnt).Data);
    
end


end


function numBins = determineBinNumber(DSet, XVars, YVars, Series, iX)
% Work out how many bins there are in the prebinned data for one participant

% Do calcuation for each series and participant. Note the code can 
% only deal with
% the case where there are the same number of bins for each
% series, and participant
numBins = NaN(length(Series), length(DSet.P));

for iPtpnt = 1 : length(DSet.P)
    
    for iY = 1 : length(YVars)
        
        xValues = XVars(iX).ProduceVar(DSet.P(iPtpnt).Data);
        
        for iS = 1 : length(Series)
            
            includedData = findIncludedTrials(DSet, iPtpnt, XVars, iX, ...
                YVars, iY, Series, iS);
            
            bins = unique(xValues(includedData));
            bins(isnan(bins)) = [];
            
            numBins(iS, iPtpnt) = length(bins);
            
        end
        
    end
    
end
                
numBins = unique(numBins);

if length(numBins) ~= 1
    error(['code can only deal with ', ...
        'the case where there are the same number of bins for each', ...
        'series.'])
end

end
   

