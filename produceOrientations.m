function orientations =  produceOrientations (nTrials, mu_cat1, mu_cat2, kappa_s, Target)

orientations = zeros(nTrials, 1);

for i = 1:nTrials
    if Target(i,1) == 0 % if category is one then compute the orientation based on cat 1 mean as follows
        orientations(i,1) =  circ_vmrnd_fixed(mu_cat1, kappa_s, 1); %generate trial orientation .
    else
        orientations(i,1) =  circ_vmrnd_fixed(mu_cat2, kappa_s, 1); %drawn from vonMises distribution.

    end
end

end


