function optimOptionsOnly = getOptimoptionsOnlyOptions(solvers)

%GETOPTIMOPTIONSONLYOPTIONS Return options that can only be set by OPTIMOPTIONS
% 
%   OPTIONS = GETOPTIMOPTIONSONLYOPTIONS returns a cell array of options
%   that can only be set via OPTIMOPTIONS.
%
%   OPTIONS = GETOPTIMOPTIONSONLYOPTIONS(SOLVERS) returns a cell array of
%   options that can only be set via OPTIMOPTIONS for the specified
%   solvers.
%
%   Private to OPTIMSET and OPTIMTOOL.

%   Copyright 2013-2016 The MathWorks, Inc.

% Determine the solvers we're going to find optimoptions-only options for
if nargin == 0
    solvers = fieldnames(createSolverStruct({'optim'}));
end

% Get all the options from all the solvers via optimoptions
numSolvers = length(solvers);
allOptimOptions = cell(1, numSolvers);
for i = 1:length(solvers)
   
    thisOpts = optimoptions(solvers{i});
    
    % Need all public properties, including hidden ones.
    mc = metaclass(thisOpts);
    numProps = length(mc.PropertyList);
    isPublicProp = false(numProps, 1);
    allNames = cell(1, numProps);
    for j = 1:numProps
        allNames{j} = mc.PropertyList(j).Name;
        isPublicProp(j) = strcmp(mc.PropertyList(j).SetAccess, 'public') && strcmp(mc.PropertyList(j).GetAccess, 'public');        
    end
    allOptimOptions{i} = allNames(isPublicProp);
    
end
allOptimOptions = unique([allOptimOptions{:}]);

% Get the options from optimset
allOptimset = fieldnames(optimset);

% Fields only supported by optimoptions
optimOptionsOnly = setdiff(allOptimOptions, allOptimset);

% Exclude options that are only in optimoptions but are deprecated and not
% yet removed
optionsDeprecatedAndNotRemoved = {};
optimOptionsOnly = setdiff(optimOptionsOnly, optionsDeprecatedAndNotRemoved);