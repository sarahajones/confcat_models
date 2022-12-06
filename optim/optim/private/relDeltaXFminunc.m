function relDelta = relDeltaXFminunc(x,deltaX)
%relDeltaXFminunc Helper function for fminunc.
% 
% Calculates the relative displacement, as used in fminunc medium-scale
% stopping test.

%   Copyright 2011 The MathWorks, Inc.

relDelta = norm(deltaX ./ (1 + abs(x)),inf);
