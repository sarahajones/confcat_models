%PROB2STRUCT Convert an OptimizationProblem to solver form. 
%
% PROBLEMSTRUCT = PROB2STRUCT(PROB) converts the OptimizationProblem input
% PROB and returns a structure that contains the solver name, solver
% options and problem description in a form the solver can use.
%
% PROBLEMSTRUCT = PROB2STRUCT(PROB, X0) additionally specifies an initial
% point. X0 is a structure with fieldnames that match the 'Name' property
% of the OptimizationVariables in the problem.
%
% See also OPTIMPROBLEM 

 
%   Copyright 2017 The MathWorks, Inc.

