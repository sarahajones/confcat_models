function isValid = validOptions(options,optionsFieldnames,minNumberOfOptimOptions)
%VALIDOPTIONS Checks if 'options' is a valid optimization option structure.
%   The input argument 'options' is to be validated
%   The input argument 'optionsFieldnames' is cell array of fields that 
%   'options' will be checked against. The third optional argument 
%   'minNumberOfOptimOptions' is the minimum number of fields to be checked.
%   The output 'isValid' is a boolean.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2012 The MathWorks, Inc.

if nargin < 3
    % By default we check for five matching fields (to make sure that
    % MATLAB optimization solvers get through this test)
    minNumberOfOptimOptions = 5; % Display, MaxFunEvals, MaxIter, TolX, and TolFun.
end
% options must be a structure
isValid = isa(options,'struct') || isa(options,'optim.options.SolverOptions');
% options field names to be validated
if isValid
    ofieldnames = fieldnames(options);
    isValid = nnz(ismember(optionsFieldnames, ofieldnames)) >= minNumberOfOptimOptions;
end
