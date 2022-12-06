function needsAt = checkFunctionField(userInput)
% CHECKFUNCTIONFIELD checks if the userInput needs a leading '@' to be used
% as a function handle

%   Private to OPTIMTOOL
%   Copyright 2014 The MathWorks, Inc.

try
  % If the userInput is a char but not a handle, we can add @.
    command0 = sprintf('ischar(%s)',userInput);
    command1 = sprintf('~isa(%s,''function_handle'')',userInput);
    needsAt = evalin( 'base', command0 ) && evalin( 'base', command1 );
catch
    needsAt = true; % default action is to add @
end

