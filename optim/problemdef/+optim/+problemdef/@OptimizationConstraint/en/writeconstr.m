%WRITECONSTR Write an OptimizationConstraint to a text file
%
%   WRITECONSTR(CON) writes a simplified mathematical form of the
%   OptimizationConstraint into a text file. The file name is the workspace
%   name of the OptimizationConstraint CON, appended with '.txt'. If
%   writeconstr cannot construct the file name from the input constraint,
%   it writes to the file 'WriteConstrOutput.txt'. writeconstr overwrites
%   any existing file.
%
%   WRITECONSTR(CON, FILENAME) writes a simplified mathematical form of the
%   OptimizationConstraint in the specified file.
%
%   See also SHOWCONSTR

 
%   Copyright 2017 The MathWorks, Inc.

