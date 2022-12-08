%script to analyse the model fits achieved by running the wrapIBSBADSFIT
%function - the output structure is loaded and analysed accordingly. 

%load in the data
data = load('ModelFit.mat');

%check all tryCounts are at 1 (they are above 1 the model fit had
%difficulty converging

bestFits = zeros(52,16); %adjust size for new bins (numPP,numModel*numParams)
averageModelFit = zeros(3,1); %adjust for number of models
for iParticipant = 1:52 %how many participants do we have
    for jModel = 1:3 %adjust for number of models
        likeli = zeros(10,1); %adjust for number of runs
        for kRun = 1:10 %adjust for number of runs
            if data.ModelFit.P(iParticipant).Model(jModel).run(kRun).try >1
                 warning('TryCount for participant %i is above 1.', iParticipant); 
            end
          likeli(kRun,1) =  data.ModelFit.P(iParticipant).Model(jModel).run(kRun).result(1,8); %store the likelihood on that run 
        end
        [val, idx] = min(likeli); %find the overall min likelihood for that model, for that participant
        bestFits(iParticipant, ((((jModel-1)*8) + 1 ):((jModel-1)*8) + 8)) = data.ModelFit.P(iParticipant).Model(jModel).run(idx).result;
        
       averageModelFit(jModel,1) = mean(bestFits(:,(jModel*8)));
       [bestFit, bestModel] = min(averageModelFit);
       
    end
end

save('bestFits.mat')

%now have the best model fit overall (model 1?)
%now within each pp, subtract their model1 fit from the others 
%put rezeroed fits into a neat matrix. 
modelFits = zeros(53,3);
for iParticipant = 1:52
    for jModel = 1:3
        modelFits(iParticipant, jModel) = bestFits(iParticipant, jModel*8) -  bestFits(iParticipant, bestModel*8);
    end
end

%calculate bootci for each model 
ci1 = bootci(10000,@mean,modelFits((1:52),1));
ci2 = bootci(10000,@mean,modelFits((1:52),2));
%add back in for further models as needed 
ci3 = bootci(10000,@mean,modelFits((1:52),3)); 
%ci4 = bootci(10000,@mean,modelFits((1:52),4));

%reaverage across the other models
for i = 1:3
modelFits(53,i) = mean(modelFits((1:52),i));
end

%plot bar of fits - are the error bars right here??
models = 1:3;
data = (modelFits(53,:))';

errlow = [abs(data(1)- ci1(1)), abs(data(2)-ci2(1)), abs(data(3)-ci3(1)) ]; % abs(data(3)-ci3(1)), abs(data(4)-ci4(1))];
errhigh = [abs(data(1)- ci1(2)), abs(data(2)-ci2(2)), abs(data(3)-ci3(2))]; % abs(data(3)-ci3(2)), abs(data(4)-ci4(2))];
errBottom = [(data(1)- ci1(1)), (data(2)-ci2(1)), (data(3)-ci3(1))]; %, (data(3)-ci3(1)), (data(4)-ci4(1))];

figure
bar(models,data)                

hold on

er = errorbar(models,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold on

colormap winter
for iPtpnt = 1:52
           plot(modelFits(iPtpnt, :), 'LineWidth', 1)
   
end

hold off


%plot out param values per participant per model
%looking at the sigma values 
for jModel = 1:2
figure
for iPtpnt = 1 : 52
           plot(1:5, bestFits(iPtpnt, ((((jModel-1)*8) + 2):((jModel-1)*8) + 6 )),'Color', [0, 0, 0], 'LineWidth', 1) %black 
           hold on
           plot(1:5, bestFits(iPtpnt, ((((jModel-1)*8) + 7):((jModel-1)*8) + 11 )), 'Color', [0.2, 0.7, 0.7], 'LineWidth', 1) %green
           hold on   
end
hold off
end


%check and see how many pp had model2 as the best fitting model
%see the variation in best fitting models across pp. 
winningModels = zeros(13,5);
for iParticipant = 1:52
    for jModel = 1:20
        winningModels(iParticipant, jModel) = bestFits(iParticipant, jModel*17);
    end
end

modelNumber = zeros(52,1);
for iParticipants = 1:52
    winningModels(iParticipants, 3) = min(winningModels(iParticipants, 1:2));
    modelNumber(iParticipants,1) = find(winningModels(iParticipants, 1:2)== winningModels(iParticipants, 3));
end

modelNumber = modelNumber';

model = zeros(4,1);
for i = 1:4
    model(i,1) = sum(modelNumber == i);
    
end

model = model';

figure
bar(models,model)