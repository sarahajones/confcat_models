function con = optimconstr(varargin)
%OPTIMCONSTR Create an optimization constraint array
%
%   CONSTR = OPTIMCONSTR(N) creates an N-by-1 vector of optimization
%   constraints.
%
%   CONSTR = OPTIMCONSTR(N1,N2,...) or OPTIMEXPR([N1 N2 ...]) creates an
%   N1-by-N2-by-... array of uninitialized optimization constraints.
%
%   CONSTR = OPTIMCONSTR(CSTR) creates an array of uninitialized
%   optimization constraints from the cell array of character vectors or
%   string vector CSTR. CONSTR is 1-by-N if CSTR is 1-by-N, otherwise
%   CONSTR is N-by-1, where N is the number of elements of CSTR. CONSTR can
%   be indexed by the character vectors or strings in CSTR.
%
%   CONSTR = OPTIMCONSTR(CSTR1,CSTR2,...) or OPTIMCONSTR({CSTR1,
%   CSTR2,...}) creates an NCSTR1-by-NCSTR2-by-... array of uninitialized 
%   optimization constraints. CSTRi is either a cell array of character
%   vectors or a string vector. NCSTRi is the number of character vectors
%   or strings in each input, CSTRi. CONSTR can be indexed in the i-th
%   dimension by the character vectors or strings in CSTRi.
%   
%   CONSTR = OPTIMCONSTR(CSTR1,N2,...) or OPTIMCONSTR({CSTR1, N2,...})
%   creates an NCSTR1-by-N2-by-... array of uninitialized optimization
%   constraints defined by a combination of numeric values and character
%   arrays or stings. CSTRi is either a cell array of character vectors or
%   a string vector. Ni is a double indicating the size of the constraints
%   in the i-th dimension.
%
%   Examples:
%   % Create a 4-by-1 constraint array
%   con1 = optimconstr(4)
%
%   % Create a 3-by-2 constraint array
%   con2 = optimconstr(3, 2)
% 
%   % Create a 1-by-3 constraint array where the columns can be 
%   % indexed by character vectors
%   con3 = optimconstr({'bos', 'jfk', 'lax'})
% 
%   % Create a 2-by-3 constraint array with character vector indexing
%   % in each direction
%   con4 = optimconstr({'corn', 'wheat'}, {'bos', 'jfk', 'lax'})
%
%   % Create a 2-by-4-by-3 array with character vector indexing in the
%   % dimensions 1 and 3
%   con5 = optimconstr({'corn', 'wheat'}, 4, {'bos', 'jfk', 'lax'})

%   See also OPTIMVAR, OPTIMEXPR, SHOWCONSTR

%   Copyright 2017 The MathWorks, Inc.

% Parse inputs
[outNames, outSize, NVpair] = optim.internal.problemdef.formatDimensionInput(varargin);

% There are no valid name value pairs for optimconstr. If anyone tries to
% use one, we need to throw a helpful error message.
if ~isempty(NVpair)
    error(message('optim_problemdef:validateIndexNames:InvalidDimensionInput'));
end

% floor size values to zero
outSize(outSize < 0) = 0;

% Forward to constraint constructor
con = optim.problemdef.OptimizationConstraint(outSize);

% Set the index names
con.IndexNames = outNames;