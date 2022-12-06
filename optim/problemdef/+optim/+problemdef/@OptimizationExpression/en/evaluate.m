%EVALUATE Compute the value of an OptimizationExpression
%
%   VALUE = EVALUATE(EXPR, VARVAL) returns the value of the
%   OptimizationExpression EXPR evaluated at the point VARVAL. VARVAL is a
%   structure with fieldnames that match the 'Name' property of the
%   OptimizationVariables composing the expression EXPR (as returned by the
%   SOLVE method of the OptimizationProblem class).
%
% Examples:
%   
%      % Create two variables and an expression to be evaluated:
%      x = optimvar('x', 1, 3);
%      y = optimvar('y', 1, 4);
%      constr = sum(x) + sum(y);
%
%      % Build a structure to represent the point at which to evaluate:
%      sol.x = [1 5 2];
%      sol.y = [0 2 1 4];
%
%      % Evaluate the expression: 
%      val = evaluate(total, sol)
%

 
%   Copyright 2017-2018 The MathWorks, Inc.

