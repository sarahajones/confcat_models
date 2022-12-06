function clearProblemFields(currentSolver)
% CLEARPROBLEMFIELDS clear problem fields in OPTIMTOOL
%   Private to OPTIMTOOL

%   Copyright 2006 The MathWorks, Inc.

% createProblemStruct first two arguments are solverName and defaultSolver.
% In this case, we want to set solverName to "currentSolver" so that all
% fields for the current solver is reset. 

if isappdata(0,'optimTool_Problem_Data')
   currentProbStruct = getappdata(0,'optimTool_Problem_Data');
end
tempstruct = createProblemStruct(currentSolver,[]);
tempstruct = rmfield(tempstruct,'solver');
fieldsToReset = fieldnames(tempstruct); 

% Reset all fieldnames 'fieldsToReset' in 'currentProbStruct' to []
for i = 1:length(fieldsToReset)
    currentProbStruct.(fieldsToReset{i}) = [];
end

setappdata(0,'optimTool_Problem_Data',currentProbStruct);

if isappdata(0,'optimTool_Problem_HashTable')
    problemHash = getappdata(0,'optimTool_Problem_HashTable');
    % Set all fields in the hash table corresponding to the problem
    % structure to empty
    for i = 1:length(fieldsToReset)
        problemHash.remove(fieldsToReset{i});
    end
end
   