%SOLVE Solve optimization problem
%
% SOL = solve(PROB) solves the specified optimization problem. The solution
% is returned as a struct.
%
% SOL = solve(PROB, X0) solves the optimization problem using the
% specified initial point. X0 is a structure with field names that match the 
% 'Name' property of the OptimizationVariables in PROB.
%
% SOL = solve(___, NAME, VALUE) specifies additional options for solve
% using one or more name-value pair arguments:
% 
%  'Solver' - Name of the solver to be used ('linprog' or 'intlinprog')
% 'Options' - Options for the solver 
% 
% [SOL,FVAL] = solve(PROB,...) returns the value of the objective at the
% solution X. FVAL matches the structure of the Objective property in PROB.
%
% [SOL,FVAL,EXITFLAG] = solve(PROB,...) returns an EXITFLAG that describes
% the exit condition. See the documentation for the possible values of
% EXITFLAG and the corresponding exit conditions.
%
% [SOL,FVAL,EXITFLAG,OUTPUT] = solve(PROB,...) returns additional
% information about the solution provided by the solver in the structure
% OUTPUT. See the documentation for OUTPUT.solver to get a complete
% description.
%
% [SOL,FVAL,EXITFLAG,OUTPUT,LAMBDA] = solve(PROB,...) returns the Lagrange
% multipliers from the solver in the LAMBDA structure. See the
% documentation for OUTPUT.solver to get a complete description.
%
% Example: 
%
%      % Create a simple optimization problem
%      x = optimvar('x', 1, 3, 'Type', 'integer', 'LowerBound', 0);
%      y = optimvar('y', 1, 4, 'LowerBound', 0);
%      p = optimproblem('ObjectiveSense', 'maximize');
%      e = sum(x) + sum(y);
%      p.Objective = e;
%      p.Constraints = [e <= 4.5; sum(y) <= 2.5];
%
%      % Solve the problem
%      val1 = solve(p)
%
%      % Build a structure to represent the initial point 
%      x0s.x = [1 1 1];
%      x0s.y = [0 0 0 0];
%
%      % Solve the problem with the initial point
%      val2 = solve(p, x0s)  
% 
% See also OPTIMPROBLEM, OPTIMOPTIONS

     
%   Copyright 2017-2018 The MathWorks, Inc.

