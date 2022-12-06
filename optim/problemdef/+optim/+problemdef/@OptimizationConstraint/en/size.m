%SIZE   Size of an OptimizationConstraint.
%   D = SIZE(CON), for M-by-N OptimizationConstraint CON, returns the
%   two-element row vector D = [M,N] containing the number of rows and
%   columns in the matrix. For N-D OptimizationConstraints, SIZE(CON)
%   returns a 1-by-N vector of dimension lengths. Trailing singleton
%   dimensions are ignored.
%
%   [M,N] = SIZE(CON) for a matrix OptimizationConstraint, CON, returns the
%   number of rows and columns in X as separate output variables.
%
%   [M1,M2,M3,...,MN] = SIZE(CON) for N > 1 returns the sizes of the first
%   N dimensions of the OptimizationConstraint CON. If the number of output
%   arguments N does not equal NDIMS(CON), then for:
%
%   N > NDIMS(CON), SIZE returns ones in the "extra" variables, 
%                 i.e., outputs NDIMS(X)+1 through N.
%   N < NDIMS(CON), MN contains the product of the sizes of dimensions N
%                 through NDIMS(CON).
%
%   M = SIZE(CON,DIM) returns the length of the dimension specified by the
%   scalar DIM.  For example, SIZE(CON,1) returns the number of rows. If
%   DIM > NDIMS(CON), M will be 1.
%
%   See also LENGTH, NDIMS, NUMEL.

 
%   Copyright 2017 The MathWorks, Inc.

