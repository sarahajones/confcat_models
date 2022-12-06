%% contrast_accuracy_trialtype_plot

%uses functions from mT folder to plot contrast level of trials (prebinned) against
%accuracy reports
%split by trial type

%runs using code by:
% Joshua Calder-Travis, j.calder.travis@gmail.com

%adapted by:
% Sarah Ashcroft-Jones 
% sarahashjones@gmail.com
% GitHub: sarahajones

%% 
XVars.ProduceVar = @(Data) Data.ContrastLevel;
XVars.ProduceVar = @(Data) Data.ContrastLevel;
XVars.NumBins = 'prebinned';

YVars.ProduceVar = @(Data, inclTrials) (mean(Data.Accuracy(inclTrials)).*100);

YVars.FindIncludedTrials = @(Data, inclTrials) true;

Series(1).FindIncludedTrials = @(Data)   Data.numGabors == 1;
Series(2).FindIncludedTrials = @(Data)   Data.numGabors == 2; 


PlotStyle.Xaxis(1).Title = 'ContrastLevel';
PlotStyle.Yaxis(1).Title = 'Percent accuracy';

PlotStyle.Data(1).Name = 'numGabors = 1';
PlotStyle.Data(1).PlotType = 'line';

PlotStyle.Data(2).Name = 'numGabors = 2';
PlotStyle.Data(2).PlotType = 'line';


figHandle = mT_plotVariableRelations(DataSet, XVars, YVars, Series, PlotStyle);