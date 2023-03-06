%script to analyse the model fits achieved by running the wrapIBSBADSFIT
%function - the output structure is loaded and analysed accordingly. 

%load in the data
data = load('ModelFit.mat');
numModel = 2;
numParticipant = 52;
%check all tryCounts are at 1 (they are above 1 the model fit had
%difficulty converging

bestFits = zeros(52,16); %adjust size for new bins (numPP,numModel*numParams)
averageModelFit = zeros(numModel,1); %adjust for number of models
for iParticipant = 1:52 %how many participants do we have
    for jModel = 1:numModel %adjust for number of models
        likeli = zeros(10,1); %adjust for number of runs
        for kRun = 1:10 %adjust for number of runs
            if data.ModelFit.P(iParticipant).Model(jModel).run(kRun).try >1
                 warning('TryCount for participant %i is above 1.', iParticipant); 
            end
          likeli(kRun,1) =  data.ModelFit.P(iParticipant).Model(jModel).run(kRun).result(1,9); %store the likelihood on that run 
        end
        [val, idx] = max(likeli); %find the overall min likelihood for that model, for that participant
        bestFits(iParticipant, ((((jModel-1)*9) + 1 ):((jModel-1)*9) + 9)) = data.ModelFit.P(iParticipant).Model(jModel).run(idx).result;
        
       averageModelFit(jModel,1) = mean(bestFits(:,(jModel*9)));
       [bestFit, bestModel] = max(averageModelFit);
       
    end
end

save('bestFits.mat')

%now have the best model fit overall (model 1?)
%now within each pp, subtract their model1 fit from the others 
%put rezeroed fits into a neat matrix. 
modelFits = zeros(53,numModel);
for iParticipant = 1:52
    for jModel = 1:numModel
        modelFits(iParticipant, jModel) = bestFits(iParticipant, bestModel*9) - bestFits(iParticipant, jModel*9);
    end
end

%calculate bootci for each model 
ci1 = bootci(10000,@mean,modelFits((1:52),1));
ci2 = bootci(10000,@mean,modelFits((1:52),2));
%add back in for further models as needed 
%ci3 = bootci(10000,@mean,modelFits((1:52),3)); 
%ci4 = bootci(10000,@mean,modelFits((1:52),4));

%reaverage across the other models
for i = 1:numModel
modelFits(53,i) = mean(modelFits((1:52),i));
end

%plot bar of fits - are the error bars right here??
models = 1:numModel;
data = (modelFits(53,:))';

errlow = [abs(data(1)- ci1(1)), abs(data(2)-ci2(1))]; %, abs(data(3)-ci3(1)) ]; % abs(data(3)-ci3(1)), abs(data(4)-ci4(1))];
errhigh = [abs(data(1)- ci1(2)), abs(data(2)-ci2(2))];% , abs(data(3)-ci3(2))]; % abs(data(3)-ci3(2)), abs(data(4)-ci4(2))];
%errBottom = [(data(1)- ci1(1)), (data(2)-ci2(1))]; %, (data(3)-ci3(1))]; %, (data(3)-ci3(1)), (data(4)-ci4(1))];

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


%compute BIC from the log likelihoodd 
%BIC: ?2 × log-likelihood + log(Number of observations) × number of parameters
numParams = 17;
numObs = 240*numParticipant;

modelBIC = zeros(numParticipant,numModel);
for iParticipant = 1:numParticipant
    for jModel = 1:numModel
        modelBIC(iParticipant, jModel) = ( -2*(bestFits(iParticipant,jModel*9))) + (log(numObs)*numParams);
    end
end


avgbic = zeros(1,numModel);
for jModel = 1:numModel
  avgbic(1, jModel) = ( -2*(averageModelFit(jModel,1))) + (log(numObs)*numParams);
end

[~,bic] = aicbic(averageModelFit,numParams,numObs);

%now have the best model fit overall (model 2)
%now within each pp, subtract their model2 fit from the others 
%put rezeroed fits into a neat matrix. 
deltaBIC = zeros(numParticipant+1,numModel);
for iParticipant = 1:numParticipant
    for jModel = 1:numModel
        deltaBIC(iParticipant, jModel) =   modelBIC(iParticipant, jModel) - modelBIC(iParticipant, bestModel);
    end
end

%deltaBIC = -1.*deltaBIC; %reverse the yaxis

%calculate bootci for each model 
ci1 = bootci(10000,@mean,deltaBIC((1:numParticipant),1));
ci2 = bootci(10000,@mean,deltaBIC((1:numParticipant),2));
%ci3 = bootci(10000,@mean,deltaBIC((1:13),3));
%ci4 = bootci(10000,@mean,deltaBIC((1:13),4));

%reaverage across the other models
for i = 1:numModel
deltaBIC(numParticipant+1,i) = mean(deltaBIC((1:numParticipant),i));
end


%plot bar of fits - are the error bars right here??
models = 1:numModel;
data = (deltaBIC(numParticipant+1,:))';

errlow = [abs(data(1)- ci1(1)), abs(data(2)-ci2(1))]; %, abs(data(3)-ci3(1)), abs(data(4)-ci4(1))];
errhigh = [abs(data(1)- ci1(2)), abs(data(2)-ci2(2))]; %, abs(data(3)-ci3(2)), abs(data(4)-ci4(2))];
errBottom = [(data(1)- ci1(1)), (data(2)-ci2(1))]; %, (data(3)-ci3(1)), (data(4)-ci4(1))];

figure
bar(models,data)                

hold on

er = errorbar(models,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold on

colormap winter
for iPtpnt = 1:numParticipant
           plot(deltaBIC(iPtpnt, :), 'LineWidth', 1)
   
end

hold off


saveas(gcf, 'BIC_fit_plot.png')


% %plot out param values per participant per model
% %looking at the sigma values 
% for jModel = 1:2
% figure
% for iPtpnt = 1 : 52
%            plot(1:5, bestFits(iPtpnt, ((((jModel-1)*8) + 2):((jModel-1)*8) + 6 )),'Color', [0, 0, 0], 'LineWidth', 1) %black 
%            hold on
%            plot(1:5, bestFits(iPtpnt, ((((jModel-1)*8) + 7):((jModel-1)*8) + 11 )), 'Color', [0.2, 0.7, 0.7], 'LineWidth', 1) %green
%            hold on   
% end
% hold off
% end


%check and see how many pp had model2 as the best fitting model
%see the variation in best fitting models across pp. 
winningModels = zeros(52,numModel+1);
for iParticipant = 1:52
    for jModel = 1:numModel
        winningModels(iParticipant, jModel) = bestFits(iParticipant, jModel*9);
    end
end

modelNumber = zeros(52,1);
for iParticipants = 1:52
    winningModels(iParticipants, numModel+1) = max(winningModels(iParticipants, 1:2));
    modelNumber(iParticipants,1) = find(winningModels(iParticipants, 1:2)== winningModels(iParticipants, numModel+1));
end

modelNumber = modelNumber';

model = zeros(numModel,1);
for i = 1:numModel
    model(i,1) = sum(modelNumber == i);
    
end

model = model';

figure
bar(models,model)