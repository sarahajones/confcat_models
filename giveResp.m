function resp = giveResp(nTrials, Data, midline, prior, lapseRate)
if prior ~= 0.5; error ('Prior is not neutral cannot compute'); end

% What response is given in each case? resp = 0-> cat 1, resp = 1 -> cat 2.
resp = zeros(nTrials, 1); %default on zero, if cat 1 percept it stays as the default unchanged

lapse = (rand(1,nTrials))';
vector1 = lapse >= lapseRate; %these trials should be computed properly
vector2 = lapse < lapseRate; %theses trials should be lapses

respIfNotLapse = nan(nTrials, 1);
respIfNotLapse(Data.Percept > (midline)) = 1;
respIfNotLapse(Data.Percept < (midline)) = 0;
respIfNotLapse(Data.Percept == (midline)) = randi([0 1], 1, 1);

if (sum(Data.Percept== (midline))) >1 
    warning ('More than 1 percept value on the category bound, check!'); 
end

assert(~any(isnan(respIfNotLapse)))
resp(vector1) = respIfNotLapse(vector1);

resp(vector2, 1) = randi([0 1], 1, sum(vector2)); %lapse out on the other trials (randomly pick or choose)

end