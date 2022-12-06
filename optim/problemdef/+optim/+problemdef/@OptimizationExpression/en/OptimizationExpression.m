classdef OptimizationExpression< matlab.mixin.CustomDisplay
%OPTIMIZATIONEXPRESSION Expressions for optimization problems
%
%   OptimizationExpression contains one or more mathematical expressions
%   composed of OptimizationVariables.
%
%   Construct an objective expression by performing mathematical operations
%   on optimization variables. For example,
%
%   x = optimvar('x');
%   y = optimvar('y');
%   expr = 2*x + 3*y
%
%   Use the OPTIMEXPR function to create an array of all zero
%   OptimizationExpressions.
%
%   See also OPTIMEXPR

     
    %   Copyright 2017-2018 The MathWorks, Inc.

    methods
        function ctranspose(in) %#ok<MANU>
            %' Transpose OptimizationExpression.
            %
            %   EXPR' transposes the OptimizationExpression, EXPR.
        end

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

        function horzcat(in) %#ok<MANU>
            %HORZCAT Horizontal concatenation of OptimizationExpressions.
            %
            %   See also HORZCAT.
        end

        function isempty(in) %#ok<MANU>
            %ISEMPTY True for empty OptimizationExpression array.
            %
            %   See also ISEMPTY.
        end

        function isscalar(in) %#ok<MANU>
            %ISSCALAR True if OptimizationExpression is scalar.
            %
            %   See also ISSCALAR.
        end

        function ldivide(in) %#ok<MANU>
            %.\ LDIVIDE Left array divide for OptimizationExpressions.
            %
            %   See also .\, LDIVIDE.
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

        function length(in) %#ok<MANU>
            % LENGTH(EXPR) Length of OptimizationExpression.
            % array.
            %
            %   See also LENGTH.
        end

        function lt(in) %#ok<MANU>
            %< LT Less than for OptimizationExpressions.
            %
            %   E1 < E2 is not supported for OptimizationExpressions.
            %
            %   Use <= instead.
        end

        function minus(in) %#ok<MANU>
            %- MINUS Subtraction of OptimizationExpressions.
            %
            %   See also -, MINUS.
        end

        function mldivide(in) %#ok<MANU>
            %\  MLDIVIDE Backslash or left matrix divide for OptimizationExpressions.
            %
            %   Backslash is used to create EOUT = X\EIN, where X must be a
            %   numeric scalar and EOUT,EIN are OptimizationExpressions.
            %
            %   See also \, MLDIVIDE.
        end

        function mrdivide(in) %#ok<MANU>
            %/  MRDIVIDE Slash or right matrix divide for OptimizationExpressions.
            %
            %   Slash is used to create EOUT = EIN/Y, where Y must be a
            %   numeric scalar and EOUT,EIN are OptimizationExpressions.
            %
            %   See also /, MRDIVIDE.
        end

        function mtimes(in) %#ok<MANU>
            %*  MTIMES Matrix multiply for OptimizationExpressions.
            %
            %   See also *, MTIMES.
        end

        function numel(in) %#ok<MANU>
            %NUMEL Number of elements in an OptimizationExpression array.
            %
            %   See also NUMEL.
        end

        function plus(in) %#ok<MANU>
            %+ PLUS Addition of OptimizationExpressions.
            %
            %   See also +, PLUS.
        end

        function rdivide(in) %#ok<MANU>
            %./ RDIVIDE Right array divide for OptimizationExpressions.
            %
            %   See also ./, RDIVIDE.
        end

        function sum(in) %#ok<MANU>
            %SUM Sum of elements of OptimizationExpression array.
            %
            %   See also SUM.
        end

        function times(in) %#ok<MANU>
            %.* TIMES Array multiply for OptimizationExpressions.
            %
            %   See also .*, TIMES.
        end

        function transpose(in) %#ok<MANU>
            % .' Transpose OptimizationExpression.
            %
            %   EXPR.' transposes the OptimizationExpression, EXPR.
        end

        function uminus(in) %#ok<MANU>
            %- UMINUS Unary minus for OptimizationExpressions.
            %
            %   See also UMINUS.
        end

        function uplus(in) %#ok<MANU>
            %+ UPLUS Unary plus for OptimizationExpressions.
            %
            %   See also UPLUS.
        end

        function vertcat(in) %#ok<MANU>
            %VERTCAT Vertical concatenation of OptimizationExpressions.
            %
            %   See also VERTCAT.
        end

    end
    methods (Abstract)
    end
    properties
        IndexNames;

    end
end
