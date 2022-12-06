classdef OptimizationConstraint< matlab.mixin.CustomDisplay
%OPTIMIZATIONCONSTRAINT Constraints for optimization problems
%
%   OptimizationConstraint contains one or more constraints expressions in
%   terms of OptimizationVariables.
%
%   Construct a constraint by constraining expressions, for example
%
%   x = optimvar('x');
%   y = optimvar('y');
%   c = 2*x + 3*y <= 5
%
%   % Create 4 equality constraints
%   Z = optimvar('Z',4,3);
%   constrsum = sum(Z,2) == 1;
%
%   Use the OPTIMCONSTR function to create an uninitialized array of
%   OptimizationConstraints.
%
%   See also OPTIMCONSTR

 
%   Copyright 2017-2018 The MathWorks, Inc.

    methods
    end
    methods (Abstract)
    end
    properties
        IndexNames;

    end
end
