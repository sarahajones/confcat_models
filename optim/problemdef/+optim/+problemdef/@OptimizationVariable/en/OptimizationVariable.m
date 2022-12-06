classdef OptimizationVariable< optim.problemdef.OptimizationExpression & handle
%OPTIMIZATIONVARIABLE Variables for optimization problems
%
%   OptimizationVariable is an array of one or more parameters to be
%   determined while solving an optimization problem.
%
%   Use the OPTIMVAR function to initialize an array of
%   OptimizationVariables.
%
%   For example, to define two distinct sets of variables for an
%   optimization problem:
%
%   % Define a 2x3 binary integer variable set.
%   x = optimvar('x',2,3,'LowerBound',0,'UpperBound',1,'Type','integer');
%
%   % Define a 2x3 variable set with named indices.
%   y = optimvar('y',{'SFO','DTW'},{'BOS','LAX','JFK'},'LowerBound',ones(2,3));
%
%   See also OPTIMVAR

     
%   Copyright 2016-2018 The MathWorks, Inc.

    methods
        function eq(in) %#ok<MANU>
            %== EQ Create an OptimizationConstraint.
            %
            %   C = E1 == E2 creates an equality OptimizationConstraint.
            %
            %   Each of E1,E2 is either an OptimizationExpression,
            %   OptimizationVariable or numeric array.
            %
            %   Example:
            %   x = optimvar('x', 2, 1);
            %   c = 2*x == 3
        end

        function ge(in) %#ok<MANU>
            %>= GE Create an OptimizationConstraint.
            %
            %   C = E1 >= E2 creates an inequality OptimizationConstraint.
            %
            %   Each of E1,E2 is either an OptimizationExpression,
            %   OptimizationVariable or numeric array.
            %
            %   Example:
            %   x = optimvar('x', 2, 1);
            %   c = 2*x >= 3
        end

        function gt(in) %#ok<MANU>
            %> GT Greater than for OptimizationExpressions.
            %
            %   E1 > E2 is not supported for OptimizationExpressions.
            %
            %   Use >= instead.
        end

        function le(in) %#ok<MANU>
            %<= LE Create an OptimizationConstraint.
            %
            %   C = E1 <= E2 creates an inequality OptimizationConstraint.
            %
            %   Each of E1,E2 is either an OptimizationExpression,
            %   OptimizationVariable or numeric array.
            %
            %   Example:
            %   x = optimvar('x', 2, 1);
            %   c = 2*x <= 3
        end

        function lt(in) %#ok<MANU>
            %< LT Less than for OptimizationExpressions.
            %
            %   E1 < E2 is not supported for OptimizationExpressions.
            %
            %   Use <= instead.
        end

    end
    methods (Abstract)
    end
    properties
        LowerBound;

        Name;

        Type;

        UpperBound;

    end
end
