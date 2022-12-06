function expr = optimexpr(varargin)
%OPTIMEXPR Create an optimization expression array
%
%   EXPR = OPTIMEXPR(N) creates an N-by-1 vector of all zero optimization
%   expressions.
%
%   EXPR = OPTIMEXPR(N1,N2,...) or OPTIMEXPR([N1 N2 ...]) creates an
%   N1-by-N2-by-... array of all zero optimization expressions.
%
%   EXPR = OPTIMEXPR(CSTR) creates a N-by-1 or 1-by-N vector of all zero 
%   optimization expressions, where N is the number of elements of CSTR. EXPR
%   is 1-by-N if CSTR is 1-by-N, otherwise EXPR X is N-by-1. EXPR can be
%   indexed by the character arrays or strings in CSTR.
%
%   EXPR = OPTIMEXPR(CSTR1,CSTR2,...) OR OPTIMEXPR({CSTR1, CSTR2,...})
%   creates an NCSTR1-by-NCSTR2-by-... array of all zero optimization
%   expressions, where NCSTR1 is the number of strings in CSTR1 and so on.
%   EXPR can be indexed by the character arrays or strings in CSTR1, CSTR2,
%   ... in the corresponding dimension.
%   
%   Example: 
%
%   % Create a vector of optimization expressions indexed by a string
%   EAST = {'jfk', 'bos'};
%   WEST = {'sfo', 'lax', 'sea'};
%   UNITCOST = [1 2 3;4 5 6];
%   amt = optimvar('amt', {EAST, WEST}, 'LowerBound', 0, 'UpperBound', 100);
%   shipcost = optimexpr(EAST);
%   for n = 1:length(EAST)
%       shipcost(n) = sum(UNITCOST(n, :).*amt(n, :));
%   end
%
%   See also OPTIMVAR, SHOWEXPR

%   Copyright 2017-2018 The MathWorks, Inc.

% Find the index of the first name-value pair input

[outNames, outSize, NVpair] = optim.internal.problemdef.formatDimensionInput(varargin);

% There are no valid name value pairs for optimexpr.  If anyone tries to
% use one, we need to throw a helpful error message.
if ~isempty(NVpair)
    error(message('optim_problemdef:validateIndexNames:InvalidDimensionInput'));
end

% floor size values to zero
outSize(outSize < 0) = 0;
% Forward to expression constructor
expr = optim.problemdef.OptimizationExpression(outSize, outNames);
