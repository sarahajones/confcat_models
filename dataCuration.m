fullData = readtable("dataframe2.csv"); %load in data (clean confidence trials only)
PID = (1:height(unique(fullData(:,44)))); %figure out how many pp we have
pp = table2array(unique(fullData(:,44))); %get their actual PIDs

%what columns do we need
data(:,1) = fullData(:,44);%source PID 1:52
data(:,2) = fullData(:,36);%source confidence 1:100
data(:,3) = fullData(:,46);%source binnedConfidence 1:5
data(:,4) = fullData(:,22);%source response (1=Zap, 0=Retrieve)
data(:,5) = fullData(:,24);%source accuracy (1=correct, 0=incorrect)
data(:,6) = fullData(:,48);%source stimulus location (recentred flipped location) (so boundary is always 350)
data(:,7) = fullData(:,49);%source distanceBin 1:6
data(:,8) = fullData(:,52);%source stimulusType 1 = bad spaceship, 0=good spaceship
data(:,9) = fullData(:,53);%source distributionType 1 = narrow, 0 = wide

data = table2array(data); %convert the table to an array

%now convert the array into a datastruct
for i = 1:length(PID)
    j = pp(i,1);
   vector1 = data(:,1) == j;
   df.P(i).data.confidence = data(vector1,2);
   df.P(i).data.binnedConf = data(vector1,3);
   df.P(i).data.response = data(vector1,4);
   df.P(i).data.accuracy = data(vector1,5);
   df.P(i).data.location = data(vector1,6);
   df.P(i).data.distanceBin = data(vector1,7);
   df.P(i).data.stimulusType = data(vector1,8);
   df.P(i).data.distribution = data(vector1,9);
end

