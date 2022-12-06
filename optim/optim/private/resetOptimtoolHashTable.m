function resetOptimtoolHashTable(oldmodelName)
%RESETOPTIMTOOLHASHTABLE reset Optimtool Java hash table stored in MATLAB
%  Private to OPTIMTOOL

%   Copyright 2006 The MathWorks, Inc.

% If oldmodel is empty then return else delete it
if isempty(getappdata(0,oldmodelName))
    return;
else
    rmappdata(0,oldmodelName);
end

