function conf = computeConfidence(nTrials,Data, sigma_s1, sigma_s2, mu_cat1, mu_cat2, Decision, Model)

numerator = NaN(nTrials, 1); %set vector for numerator
denominator = NaN(nTrials,1); %set vector for denominator
posteriorRatio = NaN(nTrials, 1); %set posterior ratio model for confidence
conf = zeros(nTrials,1); %set confidence vector
metaNoise = (randn(1, nTrials)*(Data.metacognitiveNoise))'; %computes norm dist of noise and changes the std by the metacog noise

vector1 = Model ==1 & Decision ==1; 
vector2 = Model ==1 & Decision ==0;
vector3 = Model ==2 & Decision ==1; 
vector4 = Model ==2 & Decision ==0;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL 1: compute the Bayesian posteriors
%we are building a BAyesian approximation of Confidence when the 2 means
%and 2 sigmas are different (for each category)
%first if decision == 1
numerator(vector1,1) = ((1)/(sqrt(2*pi*(sigma_s1^2)))).*exp(-((Data.Percept(vector1) + Data.Location(vector1)).^2)./(2*(sigma_s1)^2));
denominator(vector1,1) = ((1)/(sqrt(2*pi*(sigma_s2^2))))*exp(-((Data.Percept(vector1) + Data.Location(vector1)).^2)./(2*(sigma_s2)^2));
posteriorRatio(vector1,1) = log(numerator(vector1,1)./denominator(vector1,1));
posteriorRatio(vector1,1) = posteriorRatio(vector1,1) + metaNoise(vector1,1);
conf(vector1,1) = (1 ./(1 + exp(posteriorRatio(vector1,1))));
  
%now if decision ==2
numerator(vector2,1) = ((1)/(sqrt(2*pi*(sigma_s2^2)))).*exp(-((Data.Percept(vector2) + Data.Location(vector2)).^2)./(2*(sigma_s2)^2));
denominator(vector2,1) = ((1)/(sqrt(2*pi*(sigma_s1^2)))).*exp(-((Data.Percept(vector2) + Data.Location(vector2)).^2)./(2*(sigma_s1)^2));
posteriorRatio(vector2,1) = log(numerator(vector2,1)./denominator(vector2,1));
posteriorRatio(vector2,1) = posteriorRatio(vector2,1) + metaNoise(vector2,1);
conf(vector2,1) = (1 ./(1 + exp(posteriorRatio(vector2,1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL 2: compute the Bayesian posteriors (SYMMETRICAL)
%we are building a BAyesian approximation of Confidence when the 2 means
%and 2 sigmas are SAME (AVERAGE OF EACH) (for each category)
avg_sigma = (sigma_s1+ sigma_s2)/2;

%first if decision == 1
numerator(vector3,1) = ((1)/(sqrt(2*pi*(avg_sigma^2)))).*exp(-((Data.Percept(vector3) + Data.Location(vector3)).^2)./(2*(avg_sigma)^2));
denominator(vector3,1) = ((1)/(sqrt(2*pi*(avg_sigma^2))))*exp(-((Data.Percept(vector3) + Data.Location(vector3)).^2)./(2*(avg_sigma)^2));
posteriorRatio(vector3,1) = log(numerator(vector3,1)./denominator(vector3,1));
posteriorRatio(vector3,1) = posteriorRatio(vector3,1) + metaNoise(vector3,1);
conf(vector3,1) = (1 ./(1 + exp(posteriorRatio(vector3,1))));
  
%now if decision ==2
numerator(vector4,1) = ((1)/(sqrt(2*pi*(avg_sigma^2)))).*exp(-((Data.Percept(vector4) + Data.Location(vector4)).^2)./(2*(avg_sigma)^2));
denominator(vector4,1) = ((1)/(sqrt(2*pi*(avg_sigma^2)))).*exp(-((Data.Percept(vector4) + Data.Location(vector4)).^2)./(2*(avg_sigma)^2));
posteriorRatio(vector4,1) = log(numerator(vector4,1)./denominator(vector4,1));
posteriorRatio(vector4,1) = posteriorRatio(vector4,1) + metaNoise(vector4,1);
conf(vector4,1) = (1 ./(1 + exp(posteriorRatio(vector4,1))));




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %MODEL 2: compute the distance to decision bound model (at balance point of
% %350) aand translate to confidence
% distance(vector3,1) = abs(Data.Percept(vector3,1) - 350);
% %need to find a way to make this distance into confidence between 0 and 1
% %distance is to the midline at 350 , the space is 0-800 so 800 untis widde,
% %the largest distance we can have is 450. 
% conf(vector3,1) = distance(vector3,1)./350 ;
% conf(vector3,1) = (conf(vector3,1).*100) + metaNoise(vector3,1);
% if conf(vector3,1) < 0
%     conf(vector3,1) = 0;
% end
% conf(vector3,1) = conf(vector3,1)./100;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %MODEL 3: compute distance to midline between means model and translate to
% %confidence
% 
% midline = abs(mu_cat1-mu_cat2);
% midline = mu_cat1 + (midline/2);
% distance(vector4,1) = abs(Data.Percept(vector4,1) - midline);
% conf(vector4,1) = distance(vector4,1)./midline ;
% conf(vector4,1) = (conf(vector4,1).*100) + metaNoise(vector4,1);
% if conf(vector4,1) < 0
%     conf(vector4,1) = 0;
% end
% 
% 
% conf(vector4,1) = conf(vector4,1)./100;

