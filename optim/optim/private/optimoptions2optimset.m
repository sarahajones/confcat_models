function optionsStruct = optimoptions2optimset(optionsObj)
%   Create a valid optimset type structure from a optimoptions object.

%   Private to OPTIMTOOL

%   Copyright 2012-2016 The MathWorks, Inc.

% If it's a struct or it's empty
if ~isa(optionsObj,'optim.options.SolverOptions')
  optionsStruct = optionsObj;
  return;
end

% Copy fields from the object to the flat structure.

% Get fields that are different from their default value.
optionsStruct = getSetByUserOptions(optionsObj);

% Map the options to those expected by the structure in the solvers and
% optimtool.
optionsStruct = mapOptionsForSolver(optionsObj, optionsStruct);

% The options structure inside optimtool does not use the Algorithm
% option for linprog, but it is used by optimset. As such, the Algorithm
% option is mapped to LargeScale here.
if strcmpi(optionsObj.SolverName, 'linprog')
    if strcmpi (optionsObj.Algorithm, 'interior-point-legacy')
        optionsStruct.LargeScale = 'on';
    end
    if isfield(optionsStruct,'Algorithm')
        optionsStruct = rmfield(optionsStruct, 'Algorithm');
    end
end

% The options structure inside optimtool does not support
% the 'sqp-legacy' algorithm. We warn and switch to 'sqp'.
if strcmpi(optionsObj.SolverName, 'fmincon') && strcmpi(optionsObj.Algorithm, 'sqp-legacy')
    helpdlg(['''sqp-legacy'' is not available in optimtool.' char(10)...
        'Setting Algorithm to ''sqp''.'],...
        'Algorithm unavailable');
    optionsObj.Algorithm = 'sqp';
    optionsStruct.Algorithm = 'sqp';
end

function optionsStruct = getSetByUserOptions(obj)
% Go into the object and retrieve the options set by the user in a
% slimmed-down struct.
optionsStruct = getSetByUserOptionsStruct(obj);

% NOTE: If TolFunValue is set, it will cause an error later in optimtool.
% There a "valid" options structure requires a struct with at least 1
% fieldname that matches any from the set created by
% "createOptionsStruct('all')". If TolFunValue is the only one, then it
% will break.
%
% For simplicity, we will _always_ ignore TolFunValue for Optim Toolbox
% solvers. The GUI does not use TolFunValue internally, and will not export
% it, so it is easiest and least confusing to simply ignore it during import.
%
% Global Optim is different in that its options objects only contain
% TolFunValue. Although the option objects store TolFunValue, the Global
% solvers _read_ TolFun. Therefore, the GUI will work with a
% "FunctionTolerance", as long as we translate the structure field name to
% be TolFun. That is done below.

% Check to see if TolFun (OptimalityTolerance) is set.
if ~isempty(optionsStruct) && ...
   any(strcmp(obj.SolverName, {'ga','gamultiobj','patternsearch','simulannealbnd'}) ) && ...
   (isfield(optionsStruct,'TolFunValue') && ~isempty(optionsStruct.TolFunValue))
    % Only FunctionTolerance is set. Write its value into a "TolFun" field
    optionsStruct.TolFun = optionsStruct.TolFunValue;
end