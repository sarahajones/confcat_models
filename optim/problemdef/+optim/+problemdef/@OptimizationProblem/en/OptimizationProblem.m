classdef OptimizationProblem< matlab.mixin.CustomDisplay & matlab.mixin.internal.Scalar
%OPTIMIZATIONPROBLEM Optimization problem specification
%
%   OptimizationProblem contains OptimizationExpressions and
%   OptimizationConstraints that constitute an optimization problem.
%
%   Use the OPTIMPROBLEM function to initialize an OptimizationProblem.
%
%   Construct an OptimizationProblem with a custom description, for example
%
%   problem = optimproblem('Description','My Example Problem');
%
%   See also OPTIMPROBLEM

     
%   Copyright 2017 The MathWorks, Inc.

    methods
    end
    methods (Abstract)
    end
    properties
        Constraints;

        Description;

        Objective;

        ObjectiveSense;

        Variables;

    end
end
