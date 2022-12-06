function options = optimset4optimtool(solver, varargin)
%

%OPTIMSET4OPTIMTOOL Create/alter optimization OPTIONS structure for
%                   OPTIMTOOL
%
%   OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, 'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   creates an optimization options structure OPTIONS in which the named
%   parameters have the specified values.  Any unspecified parameters are
%   set to [] (parameters with value [] indicate to use the default value
%   for that parameter when OPTIONS is passed to the optimization
%   function). It is sufficient to type only the leading characters that
%   uniquely identify the parameter.  Case is ignored for parameter names.
%   NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, OLDOPTS,'PARAM1',VALUE1,...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, OLDOPTS,NEWOPTS) combines an existing
%   options structure OLDOPTS with a new options structure NEWOPTS.  Any
%   parameters in NEWOPTS with non-empty values overwrite the corresponding
%   old parameters in OLDOPTS.
%
%   OPTIONS = OPTIMSET4OPTIMTOOL creates an options structure OPTIONS where
%   all the fields are set to [].
%
%   OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER) creates an options structure with
%   all the parameter names and default values relevant to the optimization
%   function named in SOLVER. For example,
%           optimset('fminbnd')
%   or
%           optimset(@fminbnd)
%   returns an options structure containing all the parameter names and
%   default values relevant to the function 'fminbnd'.
%
%   Private function for OPTIMTOOL

%   Copyright 2013-2015 The MathWorks, Inc.

% Number of optional inputs
numVarargin = length(varargin);

% Switch on each of the optimset APIs
if (nargin==1) && (ischar(solver) || isa(solver,'function_handle') )

    % OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER)

    % Create full options structure
    options = createOptionsStruct;

    % Create options object for the specified solver
    optionsObj = optimoptions(solver);
    solverOpts = fieldnames(optionsObj);
    for i = 1:length(solverOpts)
       options.(solverOpts{i}) = optionsObj.(solverOpts{i});
    end

elseif numVarargin > 0 && ischar(varargin{1})

    % OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, 'PARAM1',VALUE1,'PARAM2',VALUE2,...)

    % Find any options that are not supported by optimset but are supported
    % by optimoptions.
    [optimsetInputs, optimoptionsNames, optimoptionsValues] = ...
        analyseNameValuePairs(varargin{:});

    % Call optimset on the optimset options
    options = optimset(optimsetInputs{:});

    % Call optimoptions on the optimoptions only options.
    options = internalOptimoptions(options, solver, ...
        optimoptionsNames, optimoptionsValues);

elseif (numVarargin == 1 && isstruct(varargin{1})) || ...
       (numVarargin > 1 && isstruct(varargin{1}) && ischar(varargin{2}))

    % OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, OLDOPTS,'PARAM1',VALUE1,...)

    % If no name-value pairs are passed, then the structure is allowed to
    % contain invalid option names which will be ignored.
    STRUCT_ONLY = length(varargin) == 1;

    % Separate the old options into optimset and optimoptions-only options
    [oldopts, oldOptimoptionsNames, oldOptimoptionsValues] = ...
        separateInputs(varargin{1});

    % Call an internal version optimset on the old optimset inputs. This
    % does the following:
    % 1. Sets options to values that are supported by optimset
    % 2. Extracts optimset options that have been set to values that
    % aren't supported by optimset. We check below whether these are
    % supported by optimoptions.
    [options, oldOptsNewValNames, oldOptsNewValValues] = ...
        internalOptimset(oldopts);

    % Find any options that are not supported by optimset but are supported
    % by optimoptions in the name-value pairs
    [optimsetInputs, newOptimoptionsNames, newOptimoptionsValues] = ...
        analyseNameValuePairs(varargin{2:end});

    % Call optimset on the new optimset options
    [options, nvPairsNewValNames, nvPairsNewValValues] = ...
        internalOptimset(options, optimsetInputs{:});

    % Call optimoptions on the optimoptions only options.
    optimoptionsNames = [oldOptimoptionsNames, newOptimoptionsNames, ...
        oldOptsNewValNames, nvPairsNewValNames];
    optimoptionsValues = [oldOptimoptionsValues, newOptimoptionsValues, ...
        oldOptsNewValValues, nvPairsNewValValues];
    options = internalOptimoptions(options, solver, ...
        optimoptionsNames, optimoptionsValues, STRUCT_ONLY);

elseif numVarargin == 2 && isstruct(varargin{1}) && isstruct(varargin{2})

    %OPTIONS = OPTIMSET4OPTIMTOOL(SOLVER, OLDOPTS,NEWOPTS)

    % Separate the old options into optimset and optimoptions-only options
    [oldopts, oldOptimoptionsNames, oldOptimoptionsValues] = ...
        separateInputs(varargin{1});

    % Call optimset on the old optimset inputs
    options = optimset(oldopts);

    % Separate the new options into optimset and optimoptions-only options
    [newopts, newOptimoptionsNames, newOptimoptionsValues] = ...
        separateInputs(varargin{2});

    % Call optimset on the new optimset inputs
    options = optimset(options, newopts);

    % Call optimoptions on the old optimoptions only options.
    optimoptionsNames = [oldOptimoptionsNames, newOptimoptionsNames];
    optimoptionsValues = [oldOptimoptionsValues, newOptimoptionsValues];
    options = internalOptimoptions(options, solver, ...
        optimoptionsNames, optimoptionsValues);

elseif nargin == 0

%   OPTIONS = OPTIMSET4OPTIMTOOL

    options = createOptionsStruct;

else
    error(message('optim:optimset4optimtool:UnsupportedAPI'));
end

%%%%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function options = internalOptimoptions(options, solver, ...
    optimoptionsNames, optimoptionsValues, STRUCT_ONLY)

if nargin < 5
    STRUCT_ONLY = false;
end

% Add all the optimoptions only options to the options structure
fnames = i_getOptimoptionsOnlyOptions;
for i = 1:length(fnames)
    if ~isfield(options, fnames{i})
        options.(fnames{i}) = [];
    end
end

% Set each specified optimoptions-only option
for i = 1:length(optimoptionsNames)

    % If the value is empty, we can set it immediately
    OPTIMOPTIONSET = isempty(optimoptionsValues{i});

    % Firstly, try to set the option using optimoptions for the current
    % solver
    if ~OPTIMOPTIONSET && ~isempty(solver)
        OPTIMOPTIONSET = i_canSetOption(solver, ...
            optimoptionsNames{i}, optimoptionsValues{i});
    end

    % Secondly, if the option still hasn't been set, loop through all the
    % solvers and try to set the option for a solver
    if ~OPTIMOPTIONSET
        % Only need the solvers that support optimoptions
        allSolvers = setdiff(i_solversInOptimtool, ...
            {'fzero', 'fminsearch', 'fminbnd', 'lsqnonneg'});
        for j = 1:length(allSolvers)
            try %#ok
                OPTIMOPTIONSET = i_canSetOption(allSolvers{j}, ...
                    optimoptionsNames{i}, optimoptionsValues{i});
                if OPTIMOPTIONSET
                    break
                end
            end
        end
    end

    % Finally, either set the option in the structure or error if we cannot
    % set the option
    if OPTIMOPTIONSET
        options.(optimoptionsNames{i}) = optimoptionsValues{i};
    elseif STRUCT_ONLY && ...
            ~ismember(optimoptionsNames{i}, fieldnames(optimset))
        % There are two possible cases left here:
        % 1) The option is not in optimset or optimoptions but has been
        % specified in a structure. To mimic the behavior of optimset, the
        % invalid option is ignored and no further action taken.
        % 2) The option is in optimset and the value is not supported. In
        % this case we should throw the error below.
    else
        error(message('optim:optimset4optimtool:CannotSetOption', optimoptionsNames{i}));
    end

end

function [options, oldOptsNewValNames, oldOptsNewValValues] = ...
    internalOptimset(oldopts, varargin)

% Initialize options and find the fields we need to set
if nargin == 1
    options = optimset();
    % Note that allFields must be a row vector, as oldOptsNewValNames is
    % assumed to be a row vector by the caller.
    allFields = fieldnames(oldopts);
    allFields = allFields(:)';
else
    options = oldopts;
    allFields = varargin(1:2:end);
    for i = 1:length(allFields)
        oldopts.(allFields{i}) = varargin{2*i};
    end
end

% Call optimset to set the options and values that are supported
numFields = length(allFields);
idxFieldSet = true(numFields, 1);
for i = 1:numFields
    if isempty(oldopts.(allFields{i})) % default
        continue;
    end
    try
        options = optimset(options, allFields{i}, oldopts.(allFields{i}));
    catch
        idxFieldSet(i) = false;
    end
end

% Find the values for existing options that are not supported by optimset
oldOptsNewValNames = allFields(~idxFieldSet);
numoldOptsNewVal = length(oldOptsNewValNames);
oldOptsNewValValues = cell(1, numoldOptsNewVal);
for i = 1:numoldOptsNewVal
    oldOptsNewValValues{i} = oldopts.(oldOptsNewValNames{i});
end

function options = createOptionsStruct

% Create an empty structure containing all optimset fields and all
% optimoptions-only fields
allFields = [fieldnames(optimset); i_getOptimoptionsOnlyOptions.'];
for i = 1:length(allFields)
    options.(allFields{i}) = [];
end

function [unsupportedOpts, idx] = findOptimoptionsOnlyOpts(optsSpecified)

% Find optimoptions-only options in a list of options
allUnsupportedFnames = i_getOptimoptionsOnlyOptions;
[unsupportedOpts, ~, idx] = intersect(allUnsupportedFnames, optsSpecified);

function [inputs, unsupportedOpts, unsupportedValues] = ...
    analyseNameValuePairs(varargin)

% Split name-value pairs into options and values specified
optsSpecified = varargin(1:2:end);
valuesSpecified = varargin(2:2:end);
numOptions = length(optsSpecified);

% Find optimoptions-only options
[unsupportedOpts, idxUnsupp] = findOptimoptionsOnlyOpts(optsSpecified);
unsupportedValues = valuesSpecified(idxUnsupp);

% Remove unsupported options from call to optimset
idxSupp = setdiff(1:numOptions, idxUnsupp);
optsSpecified = optsSpecified(idxSupp);
valuesSpecified = valuesSpecified(idxSupp);

% Call optimset
numSupp = length(optsSpecified);
inputs = cell(1, 2*numSupp);
inputs(1:2:end) = optsSpecified;
inputs(2:2:end) = valuesSpecified;


function solverList = i_solversInOptimtool

% List of solvers supported in optimtool
solverList = fieldnames(createSolverStruct({'matlab', 'optim'}));

function [oldopts, oldOptimoptionsNames, oldOptimoptionsValues] = ...
    separateInputs(oldopts)

% Set the old options in the structure
foldopts = fieldnames(oldopts);

% Convert structure to name-value pairs
nvPairsOld = cell(1, 2*length(foldopts));
for i = 1:length(foldopts)
    nvPairsOld{2*i-1} = foldopts{i};
    nvPairsOld{2*i} = oldopts.(foldopts{i});
end

% Find any options that are not supported by optimset but are supported
% by optimoptions.
[~, oldOptimoptionsNames, oldOptimoptionsValues] = ...
    analyseNameValuePairs(nvPairsOld{:});

% Return the optimset-only options as a structure
oldopts = rmfield(oldopts, oldOptimoptionsNames);

function allUnsupportedOpts = i_getOptimoptionsOnlyOptions

% Get list of all solvers in Optimization toolbox
solvers = fieldnames(createSolverStruct({'optim'}));

% optimtool currently does not support intlinprog, so we do not want these
% options in the optimtool options structure
solvers = setdiff(solvers, {'intlinprog'});

% Call getOptimoptionsOnlyOptions with all the Optimization toolbox
% solvers minus intlinprog
allUnsupportedOpts = getOptimoptionsOnlyOptions(solvers);

function OPTIMOPTIONSET = i_canSetOption(solver, optionName, optionValue)

try
    optimoptions(solver, optionName, optionValue);
    OPTIMOPTIONSET = true;
catch
    % We could get here if we try to set a numerical option to a
    % special string. If the string is equal to the default value,
    % it is fine to set the numerical option to the value.
    defaultObject = optimoptions(solver);
    OPTIMOPTIONSET = ischar(optionValue) && ...
        isprop(defaultObject, optionName) && ...
                    strcmp(optionValue, defaultObject.(optionName));
end