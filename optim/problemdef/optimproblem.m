function prob = optimproblem(varargin)
%OPTIMPROBLEM Create an optimization problem. 
%
%   PROB = OPTIMPROBLEM creates an optimization problem with default
%   properties.
%
%   PROB = OPTIMPROBLEM(NAME, VALUE, ...) creates an optimization problem 
%   with specified values of optional parameters:
%
%     'ObjectiveSense'   Sense of optimization: 'maximize', 
%                        'minimize' (default), 'max', or 'min'
%     'Objective'        A scalar OptimizationExpression
%     'Constraints'      An OptimizationConstraint array or a struct with 
%                        OptimizationConstraint arrays as fields
%     'Description'      Problem description: character array or string.
%   
%   The following simple example illustrates how an optimization problem is
%   created and solved
%
%   % Create an optimization problem
%   prob = optimproblem('Description', 'Simple Example');
%  
%   % Create a 2-by-1 optimization variable
%   x = optimvar('x', 2, 1, 'LowerBound', 0);
%
%   % Define the objective
%   prob.Objective = -5*x(1) - x(2);
%
%   % Define the constraints
%   prob.Constraints.MyConstraint1 = x(1) + x(2) <= 5;
%   prob.Constraints.MyConstraint2 = 2*x(1) + 0.5*x(2) <= 8;
%
%   % Solve the problem
%   sol = solve(prob);
%
%   See also OPTIMVAR, OPTIMCONSTR, OPTIMEXPR

%   Copyright 2017 The MathWorks, Inc.

% Forward to problem constructor
prob = optim.problemdef.OptimizationProblem();

for i = 1:2:length(varargin)
    % Let MATLAB handle whether the property exists, and let the set
    % methods handle the type validity checking
    if isstring(varargin{i})
        varargin{i} = char(varargin{i});
    elseif ~ischar(varargin{i})
        error(message('optim_problemdef:OptimizationProblem:PropNameNotChar'));
    end
        
    prob.(varargin{i}) = varargin{i+1};
end







