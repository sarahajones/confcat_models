%WRITEEXPR Write an OptimizationExpression to a text file
%
%   WRITEEXPR(EXPR) writes a simplified mathematical form of the
%   OptimizationExpression into a text file. The file name is the workspace
%   name of the OptimizationExpression EXPR, appended with '.txt'. If
%   writeexpr cannot construct the file name from the input expression, it
%   writes to the file 'WriteExprOutput.txt'. writeexpr overwrites any
%   existing file.
%
%   WRITEEXPR(EXPR, FILENAME) writes a simplified mathematical form of the
%   OptimizationExpression into FILENAME. 
%
%   See also SHOWEXPR

 
%   Copyright 2017-2018 The MathWorks, Inc.

