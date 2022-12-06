function samples = circ_vmrnd_fixed(mu, kappa, shape)

% There is a bug in circ_vmrnd such that for kappa_s ~ 0 it returns values on
% the interval [0 2*pi] rather than [-pi pi]. It also does not reshape the 
% output. Correct this...

%created by JOshua Calder-Travis, j.calder-travis@gmail.com
% utilised by Sarah Ashcroft-Jones, sarah.ashcroft-jones@psy.ox.ac.uk

samples = circ_vmrnd(mu, kappa, shape);

if kappa < 1e-6
    
    samples = samples - pi;
    
    
    samples = reshape(samples, shape);
    
    
end


% Defensive programming
if any(any(samples > pi)) || any(any(samples < -pi))
    
    error('Von Misses function is returning unexpected values')
    
    
end
