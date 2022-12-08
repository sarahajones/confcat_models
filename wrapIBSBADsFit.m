function ModelFit = wrapIBSBADsFit

data = load('cleanData.mat');

for iParticipant = 1:52 %1:52 %for each participant
    for jModel = 2:3 %1:3 %for each model
        for kRun = 1:10 %1:10 %for ten runs
            isComplete = 0;
            tryCount = 1;
            while isComplete == 0  
                try
                    ModelFit.P(iParticipant).Model(jModel).run(kRun).result = runIBSBADs(iParticipant, jModel, data);
                    isComplete = 1;
                 catch
                     warning('BADs failed to converge on attempt %i.  Trying again.', tryCount);  
                     isComplete = 0;
                     tryCount = tryCount + 1;
                 end  
            end
            ModelFit.P(iParticipant).Model(jModel).run(kRun).try = tryCount;
        end
    end
end
save("ModelFit")
end

