%SIZE   Size of an OptimizationExpression.
%   D = SIZE(EXPR), for M-by-N OptimizationExpression EXPR, returns the
%   two-element row vector D = [M,N] containing the number of rows and
%   columns in the matrix. For N-D OptimizationExpressions, SIZE(EXPR)
%   returns a 1-by-N vector of dimension lengths. Trailing singleton
%   dimensions are ignored.
%
%   [M,N] = SIZE(EXPR) for a matrix OptimizationExpression, EXPR, returns
%   the number of rows and columns in X as separate output variables.
%
%   [M1,M2,M3,...,MN] = SIZE(EXPR) for N > 1 returns the sizes of the first
%   N dimensions of the OptimizationExpression EXPR. If the number of
%   output arguments N does not equal NDIMS(EXPR), then for:
%
%   N > NDIMS(EXPR), SIZE returns ones in the "extra" variables, 
%                 i.e., outputs NDIMS(X)+1 through N.
%   N < NDIMS(EXPR), MN contains the product of the sizes of dimensions N
%                 through NDIMS(EXPR).
%
%   M = SIZE(EXPR,DIM) returns the length of the dimension specified by the
%   scalar DIM.  For example, SIZE(EXPR,1) returns the number of rows. If
%   DIM > NDIMS(EXPR), M will be 1.
%
%   See also LENGTH, NDIMS, NUMEL.

 
%   Copyright 2017 The MathWorks, Inc.

