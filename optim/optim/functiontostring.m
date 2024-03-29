function strfcn = functiontostring(fcn)
%

% FUNCTIONTOSTRING Convert the function to a string for printing.

%   Copyright 1990-2006 The MathWorks, Inc.

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch
        strfcn = '(name not printable)';
    end
end
