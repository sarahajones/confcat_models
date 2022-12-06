function stop = savemilpsolutions(x,optimValues,state,varargin)
%SAVEMILPSOLUTIONS save all the integer solution found by intlinprog.
%   All the unique solutions are created as variables (xIntSol, fIntSol) in
%   the base workspace.
%
%   See also intlinprog

%   Copyright 2014 The MathWorks, Inc.

% No reason to stop the solver from this function.
stop = false;

switch (state)
  case 'init'
    % Create new (or reset) variables in the base workspace  to collect
    % integer solution.
    assignin('base','xIntSol',[])
    assignin('base','fIntSol',[])
  case 'iter'
    % 'x' is not empty only when there is a new integer solution found
    % during branch and bound phase.
    if ~isempty(x)
      % Save integer x by appending new columns in the matrix xIntSol
      xInts = evalin('base','xIntSol');
      xInts = [xInts, x(:)];
      % Save fval for integer x by appending to the vector fIntSol
      fval = optimValues.fval;
      fInts = evalin('base','fIntSol');
      fInts = [fInts, fval];
      
      % Update the variables
      assignin('base','fIntSol',fInts)      
      assignin('base','xIntSol',xInts)
    end
  case 'done'
    % Nothing to do here.
   
end
