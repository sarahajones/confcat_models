function [percept] = producePercept(Data)

noise = randn(length(Data.Sigma_X), 1);
noise = noise.*(Data.Sigma_X); %compute some noise here to add to value
percept = Data.Location + noise; 

end 