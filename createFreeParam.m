function freeParam = createFreeParam 

lapseRate = 0.25; %mid-way between plausible upper and lower bounds
sigma_X = [1, 1.25, 1.5, 1.75, 2, 1.1, 1.35, 1.6, 1.85, 2.1]; %internal noise that is contrast and numGabor dependent
metacogNoise = 0.9; %midway between plausible upper and lower bounds
thresh = [.25, .5, .75]; 

freeParam  = [lapseRate, sigma_X, metacogNoise, thresh]; 

end
