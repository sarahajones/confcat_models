 function response = passSimulation (freeParams, designMatrix)

%unpack vector of free params
freeParam.lapseRate = freeParams(1);
freeParam.sigma_X = freeParams(2);
freeParam.metacogNoise = freeParams(3);
freeParam.confLapse = freeParams(4);

for i = 5:length(freeParams)
freeParam.thresh(i-4) = freeParams(i);
end
freeParam.thresh = sort(freeParam.thresh);

%unpack design matrix 
S.Stimulus = designMatrix(:, 4);
S.nTrials = length(S.Stimulus);
S.Location = designMatrix(:, 1);
S.Distribution = designMatrix(:,3);
S.Model = designMatrix(1,5);

%set fixed parameters 
fixedParam.midline = 350; %objective midline between (in all trials as we have read in the flipped and recentered data). 
fixedParam.mu_cat1 = 222; %128 distance to bound
fixedParam.mu_cat2 = 522; %170 distance to bound 
fixedParam.std_cat1 = 80; %narrow distribution std
fixedParam.std_cat2 = 140; %wide distribution std
fixedParam.prior = 0.5; %assume neutral prior for symmetry of decisions

%structure of responses and confidence
Response = runTrialSimulation(fixedParam, S, freeParam);

%TRANSFORM response into a Matrix
response(:, 1) = Response.Decision;
response(:, 2) = Response.binnedConfidence;

 end
 