function [msg,ID] = getLastWarnForOptimtool()
%

%getLastWarnForOptimtool checks lastwarn for inactive warnings before
%returning control to the output functions controlling the solver run.
%
% [msg,ID] = getLastWarnForOptimtool
% If the last warning thrown was displayed, msg and ID will contain the
% warning text and ID. Otherwise, it will be empty.

%   Copyright 2012-2013 The MathWorks, Inc.

[msg,ID] = lastwarn;

if ~isempty(msg)
    % Find the ID of the last displayed message
    warnStruct = warning('QUERY', 'last');
    % If the lastwarn isn't the same as the last displayed message, then do
    % not bother to display
    % NOTE: watch out for a 0x0 struct from querying the last warning
    if isempty(warnStruct) || ...
            ~strcmpi(warnStruct.identifier,ID)
        msg = ''; ID = '';
    end
end
