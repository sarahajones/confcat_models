function ResponseMatrix = createResponseMatrix (ParticipantNum, data)

ResponseMatrix(:, 1) = data.df.P(ParticipantNum).data.response;
ResponseMatrix(:, 2) = data.df.P(ParticipantNum).data.binnedConf;

end
