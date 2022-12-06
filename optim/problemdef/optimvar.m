function x = optimvar(name, varargin)
%OPTIMVAR Create optimization variables   
% 
%   X = OPTIMVAR(NAME) creates a scalar optimization variable. Note that the 
%   variable is symbolic and does not contain data. The variable can 
%   be used to create linear expressions for an optimization problem. 
% 
%   X = OPTIMVAR(NAME,N) creates an N-by-1 vector of optimization variables.  
%
%   X = OPTIMVAR(NAME,CSTR) creates a N-by-1 or 1-by-N vector of optimization 
%   variables, where N is the number of elements of CSTR. X is 1-by-N if CSTR 
%   is 1-by-N, otherwise X is N-by-1. X can be indexed by the character 
%   arrays or strings in CSTR. 
% 
%   X = OPTIMVAR(NAME,N1,N2,N3,...) or OPTIMVAR(NAME,[N1 N2 N3 ...]) creates 
%   an N1-by-N2-by-N3-by-... array of optimization variables. 
% 
%   X = OPTIMVAR(NAME,CSTR1,CSTR2,...) or OPTIMVAR(NAME,{CSTR1, CSTR2, ...}) 
%   creates an %   NCSTR1-by-NCSTR2-by-... array of optimization variables,  
%   where NCSTR1 is the number of strings in CSTR1 and so on. X can be  
%   indexed by the character arrays or strings in CSTR1, CSTR2, ... in the 
%   corresponding dimension. 
% 
% Examples: 
%    
%      % Create a 3-by-1 vector of optimization variables: 
%      x = optimvar('x', 3, 1); 
% 
%      % Create a 2-by-3 matrix of optimization variables: 
%      x = optimvar('x', [2 3]); 
% 
%      % Create a 2-by-3-by-4-by-2 array of optimization variables:   
%      x = optimvar('x', 2, 3, 4, 2); % or x = optimvar('x', [2 3 4 2]);  
% 
%      % Create a 1-by-1 optimization variable that can be indexed by 
%      % 'cars' 
%      x = optimvar('x', {'cars'}); 
% 
%      % Create a 2-by-1 optimization variable where the first element can 
%      % be indexed by 'red' and the second by 'green' 
%      x = optimvar('x', {'red';'green'}); 
%      % 
%      % Show that x('green') is the second element of the array 
%      x('green')  
%      x(2) 
% 
%      % Create a 3-by-2 matrix of optimization variables indexed by: 
%      % Rows: {'coffee', 'tea', 'water'} 
%      % Columns: {'winter', 'summer'} 
%      numdrinks = optimvar('numdrinks', {'coffee', 'tea', 'water'}, ... 
%            {'winter', 'summer'}); 
%      % 
%      % Show that x('tea', 'summer') is the (2, 2)-th element 
%      numdrinks('tea', 'summer')  
%      numdrinks(2, 2) 
%       
%      % Create a 3-by-3-by-3 array of optimization variables indexed by 
%      % Dimension 1: {'shirts', 'trousers', 'socks'} 
%      % Dimension 2: {'London','Washington','Paris'} 
%      % Dimension 3: {'Birmingham','NewYork','Toulouse'} 
%      flow = optimvar('flow', {'shirts', 'trousers', 'socks'}, ... 
%            {'London','Washington','Paris'}, ... 
%            {'Birmingham','NewYork','Toulouse'}); 
%    
%   X = OPTIMVAR(..., NAME, VALUE, ...) creates an array of optimization 
%   variables with specified values of optional parameters: 
% 
%      'Type', Variable type: 'continuous' (default), 'integer' 
%      'LowerBound', Variables lower bound: numeric, default -inf(size(x)) 
%      'UpperBound', Variables upper bound: numeric, default  inf(size(x)) 
%     
%   Example: 
%      % Create a 3-by-2 matrix of positive optimization variables that can 
%      % be displayed and referenced by 'numdrinks' 
%      numdrinks = optimvar('numdrinks', {'coffee', 'tea', 'water'}, 2, ... 
%      'LowerBound', 0); 
% 
%   Once you have created some optimization variables, you can create 
%   linear expressions from them, for example: 
% 
%   % Create optimization variables 
%   amt = optimvar('amt',{'pennies', 'nickels', 'dimes'}) 
%  
%   % Create a linear expression, cost 
%   cost = amt('pennies') + 5*amt('nickels') + 10*amt('dimes')
% 
%   See also OPTIMPROBLEM, OPTIMEXPR, SHOWBOUNDS

%   Copyright 2017 The MathWorks, Inc.

if (nargin == 0)
    error(message('optim_problemdef:optimvar:NotEnoughInputs'));
end

% Disallow leading or trailing whitespace
name = strip(string(name));

% name must adhere to same rules as MATLAB variable names
if ~isvarname(name)
    error(message('optim_problemdef:optimvar:NotMATLABVarName', namelengthmax));
end

% Find the index of the first name-value pair input
ValidNVpairs = ["Type","LowerBound","UpperBound"];

[outNames, outSize, NVpair] = optim.internal.problemdef.formatDimensionInput(varargin);

if any(outSize <= 0)
    error(message('optim_problemdef:optimvar:CannotCreateEmptyOptimVar'));
end

% Forward to variable constructor
x = optim.problemdef.OptimizationVariable(name, outSize, outNames);

% set all properties
for i = 1:2:numel(NVpair)
    fieldName = validatestring(NVpair{i}, ValidNVpairs);
    x.(char(fieldName)) = NVpair{i+1};
end


