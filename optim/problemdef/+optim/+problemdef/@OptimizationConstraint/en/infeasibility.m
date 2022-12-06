%INFEASIBILITY Compute the infeasibility of an OptimizationConstraint
%
%   VAL = INFEASIBILITY(CONSTR, VARVAL) returns the infeasibility of the
%   OptimizationConstraint CONSTR at the points in VARVAL. VARVAL is a
%   structure with fieldnames that match the 'Name' property of the
%   OptimizationVariables composing the constraint CONSTR (as returned by
%   the SOLVE method of the OptimizationProblem class).
%
% Examples:
%   
%      % Create two variables and a constraint to be evaluated:
%      x = optimvar('x', 1, 3);
%      y = optimvar('y', 1, 4);
%      constr = sum(x) + sum(y) <= 10;
%
%      % Build a structure to represent the point at which to evaluate:
%      sol.x = [1 5 2];
%      sol.y = [0 2 1 4];
%
%      % Compute the infeasibility:  
%      val = infeasibility(constr, sol)
%
%   See also optim.problemdef.OptimizationExpression/evaluate

 
%   Copyright 2017 The MathWorks, Inc.

