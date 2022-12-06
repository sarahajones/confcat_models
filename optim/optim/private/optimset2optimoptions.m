function [optionsObj, optionsStruct] = optimset2optimoptions(solver,optionsStruct)
%   Create a valid optimoptions object from a optimset type structure.

%   Private to OPTIMTOOL

%   Copyright 2012-2015 The MathWorks, Inc.

% These solvers still use flat structure.
solversWithStructOpts = {'fminsearch','fzero', 'fminbnd','lsqnonneg'};

% No work to be done for flat structure case.
if any(strcmp(solver,solversWithStructOpts))
  optionsObj = optionsStruct;
  return;
end
% Initialize the return argument for the requested solver.
optionsObj = optimoptions(solver);
defaultObj = optionsObj;
if ~isempty(optionsStruct) && isa(optionsStruct,'struct')
  StructFieldNames = fieldnames(optionsStruct);
  ObjectFieldNames = getOptionNames(optionsObj);  
  
  % Copy fields from the structure to the object if they are common
  for i = 1: length(StructFieldNames)
    fname = optionsStruct.(StructFieldNames{i});
    if ~isempty(fname) && any(strcmp(StructFieldNames{i},ObjectFieldNames))
        try
            optionsObj.(StructFieldNames{i}) = optionsStruct.(StructFieldNames{i});
        catch ME
            % The value in the structure may be a special string for a
            % numerical option. If this is the case that means the option
            % takes its default value and we can just move on. Otherwise,
            % we rethrow the error.
            currValue = optionsStruct.(StructFieldNames{i});
            if ~ischar(currValue) || ~strcmpi(currValue, defaultObj.(StructFieldNames{i}))
                rethrow(ME);                
            end
       end
            
    end
  end
  
  % Solver specific mappings required by optimoptions.
  [optionsObj, optionsStruct] = mapOptimsetToOptions(optionsObj, optionsStruct);
end
