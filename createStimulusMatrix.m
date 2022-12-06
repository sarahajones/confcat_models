function S = createStimulusMatrix (ModelNum, ParticipantNum, data)
%we want: location, distribution, stimulusType,
%stitched together into a design*nTrials matrix

%add location
S(:, 1) = data.df.P(ParticipantNum).data.location;

%add binnedLocation
S(:, 2) = data.df.P(ParticipantNum).data.distanceBin;

%add distribution
S(:, 3) = data.df.P(ParticipantNum).data.distribution;
    
%add stimulus type
S(:, 4) = data.df.P(ParticipantNum).data.stimulusType;
                     
%add Model Number
S(:, 5) = ModelNum;

end