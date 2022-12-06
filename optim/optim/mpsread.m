function problem = mpsread(mpsfile)
%MPSREAD reads LP and MILP optimization data from MPS formatted file.
%
% PROBLEM = mpsread(mpsfile) mpsfile is a string containing a file name.
% The file should contain MPS formatted data. On successful file read,
% PROBLEM is a structure that can be passed directly to either intlinprog
% or linprog functions.
%
%   See also INTLINPROG, LINPROG.

%   Copyright 2015 The MathWorks, Inc.

if (isstring(mpsfile) && isscalar(mpsfile))
    mpsfile = char(mpsfile);
end

% Read MPS file 
[problem.f,intcon,problem.Aineq,problem.bineq, ...
    problem.Aeq,problem.beq,problem.lb,problem.ub] = readMPSfile(mpsfile);

% Find indices of variables with integer constraints (non-zero).
problem.intcon = find(intcon ~= 0);
if ~isempty(problem.Aineq)
    % Remove unrestricted rows i.e, constraints with infinite RHS.
    % MPS file may have these rows but they are not needed to solve problems.
    unrestricted_rows = isinf(problem.bineq);
    if nnz(unrestricted_rows) > 0
        problem.Aineq = problem.Aineq(~unrestricted_rows,:);
        problem.bineq = problem.bineq(~unrestricted_rows);
    end
end

% Add solver and options fields depending on problem type.
if isempty(problem.intcon)
    problem.solver = 'linprog';
else
    problem.solver = 'intlinprog';
end
problem.options = optimoptions(problem.solver);

