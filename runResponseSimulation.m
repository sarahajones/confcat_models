function Response = runResponseSimulation(Data, fixedParam, S, freeParam)
%% DECISION make decision / give response
Response.Decision = giveResp(S.nTrials, Data, fixedParam.midline, fixedParam.prior, freeParam.lapseRate);

%% CONFIDENCE calculate confidence value for each trial 
Data.metacognitiveNoise = freeParam.metacogNoise;
Response.Confidence = computeConfidence(S.nTrials, Data,  fixedParam.std_cat1, fixedParam.std_cat2, fixedParam.mu_cat1, fixedParam.mu_cat2, Response.Decision, S.Model);

freeParam.thresh = sort(freeParam.thresh);
freeParam.thresh = [0, freeParam.thresh, 1];
Response.binnedConfidence = discretize(Response.Confidence,  freeParam.thresh);

%overwrite some confidence reports with a random confidence report
lapse = (rand(1, S.nTrials))';
vector1 = lapse < freeParam.confLapse; %theses trials should be lapses
Response.binnedConfidence(vector1) = randi(5, sum(vector1), 1); %bin changes here too !!

end
