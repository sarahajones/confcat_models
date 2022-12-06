function [optionsStruct,optionsFcn]  = createOptionsStruct(solverName,useValues)
%

%CREATEOPTIONSSTRUCT Create options structure for different solvers
%   Create options structure for 'solverName'. If defaultSolver is [] then
%   'fmincon' is assumed to be the default solver. The optional third
%   argument is used to populate the options structure 'optionsStruct' with
%   the values from 'useValues'.
%
%   Private to OPTIMTOOL

%   Copyright 2007-2015 The MathWorks, Inc.

if nargin < 2
    useValues = [];
end
% Perform a license check for optional solvers
if ~isempty(ver('globaloptim')) && license('test','GADS_Toolbox')
    enableAllSolvers = true;
else
    enableAllSolvers = false;
end

% Call appropriate options setting function for each solver
switch solverName
    case {'fmincon','fminunc','lsqnonlin','lsqcurvefit','linprog', ...
            'quadprog','fgoalattain','fminimax','fseminf','fsolve', ...
            'lsqlin'} 
        optionsStruct = optimset4optimtool(solverName, useValues);
        optionsFcn = 'optimoptions';
    case {'fminsearch','fzero','fminbnd','lsqnonneg'}
        optionsStruct = optimset4optimtool(solverName, useValues);
        optionsFcn = 'optimset';
    case {'ga','gamultiobj'}
        optionsStruct = gaoptimset(useValues);
        optionsFcn = 'optimoptions';
    case 'patternsearch'
        optionsStruct = psoptimset(useValues);
        optionsFcn = 'optimoptions';
    case 'simulannealbnd'
        optionsStruct = saoptimset(useValues);
        optionsFcn = 'optimoptions';
    case 'all'
        allfields = fieldnames(optimset4optimtool);
        if enableAllSolvers
            data = load('OPTIMTOOL_OPTIONSFIELDS.mat','globaloptimOptions');
            allfields = [allfields; data.globaloptimOptions];
        end
        optionsStruct = createEmptyStruct(allfields);
        optionsFcn = '';
end

% Copy the values from the struct 'useValues' to 'optionsStruct'.
if ~isempty(useValues)
    % Check for recently invalidated options that will be removed from the
    % options structure. Warn the user about them if they are set.
    useValues = checkDeprecatedOptions(useValues);
    
    copyfields = fieldnames(optionsStruct);
    Index = ismember(copyfields,fieldnames(useValues));
    for i = 1:length(Index)
        if Index(i)
            optionsStruct.(copyfields{i}) = useValues.(copyfields{i});
        end
    end
end
% We know that 'Display' field is in the optionsStruct. If it is empty then
% we set it to the GUI's default ('off'). 
if isempty(optionsStruct.Display)
   optionsStruct.Display = 'off'; 
end
%-----------------------------------------------------
function optionsStruct = createEmptyStruct(allfields)

optionsStruct = struct();
for i = 1:length(allfields)
    if ~isfield(allfields{i},optionsStruct)
        optionsStruct.(allfields{i}) = [];
    end
end

%-----------------------------------------------------
function userOptions = checkDeprecatedOptions(userOptions)
% 

% Check for options that were valid but are now deprecated. These fields
% will be removed from the options struct before returning from this
% function.

deprecatedOpts = {'LevenbergMarquardt';'NonlEqnAlgorithm'};

% Create a list of options that the user has set but are no longer valid.
% These are not meant to be translated since they are option names.
idxDeprecatedOptionsSet = false(1, length(deprecatedOpts));
for k = 1:length(deprecatedOpts)
    if isfield(userOptions,deprecatedOpts{k}) && ...
            ~isempty(userOptions.(deprecatedOpts{k}))
        idxDeprecatedOptionsSet(k) = true;
    end
end
userSetOpts = deprecatedOpts(idxDeprecatedOptionsSet);


if ~isempty(userSetOpts)
    % Strjoin requires a 1xN cell array of strings
    userSetOpts = strjoin(userSetOpts(:)', ', ');
    msg = message('optim:createOptionsStruct:DeprecatedOptions',userSetOpts);
    % Get a handle to the GUI
    optimtoolGUI = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI');
    if ~isempty(optimtoolGUI)
        % Throw warning in the 'Status and Results' panel
        javaMethodEDT('appendResults',optimtoolGUI,['Warning: ',getString(msg)]); 
    else
        % Something is wrong, write warning to command window
        warning(msg);
    end
end

% Add a check for UseParallel option
 if isfield(userOptions,'UseParallel') && ~isempty(userOptions.UseParallel)
   userOptions.UseParallel = validateopts_UseParallel(userOptions.UseParallel, ...
     true,true);
 end



