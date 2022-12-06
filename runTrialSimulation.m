function [Response] = runTrialSimulation(fixedParam, S, freeParam)
%% NOISE TO GAIN MEASURE X (PERCEPT)
%apply noise to the orientation to gain percept value
Data.Sigma_X = freeParam.sigma_X;

%housekeeping to pass off to percept function
Data.ContrastLevel = 1; %hold contrast as complete as no noise is layered over the trials
Data.Location = S.Location;
Data.Percept = producePercept(Data); %calculate perceived stimulus value

% If any percepts are outside the range [-pi pi] then move them back 
%Data.Percept = vS_mapBackInRange(Data.Percept, -pi, pi); 
% removed as assuming we are approximating a linear distribution here

Response = runResponseSimulation(Data, fixedParam, S, freeParam);    
end