function simulation = setUpSimulation(S)
%% CHOOSE THE MODEL 
S.Model = 1; %CHOOSE MODEL 1, 2, 3, 4

%% CREATE THETA VECTOR OF PARAMETER VALUES
%Key values and terms
%set fixed params 
fixedParam.mu_cat1 = (1/16).*(pi);
fixedParam.mu_cat2 = (-1/16).*(pi);
fixedParam.kappa_s = 7;%%kappa (concentration parameter, needs to be convereted for derivations to sigma)
fixedParam.sigma_s = sqrt(1/7);
fixedParam.contrasts = [0.1, 0.2, 0.3, 0.4, 0.8]; %external noise
fixedParam.prior = 0.5; %assume neutral prior for symmetry of decisions

simulation = runTrialSimulation(fixedParam, S, freeParam);

