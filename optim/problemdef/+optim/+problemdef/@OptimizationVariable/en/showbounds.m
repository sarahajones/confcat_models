%SHOWBOUNDS Display bounds on optimization variables
%
%   SHOWBOUNDS(X) displays all the finite bounds for each variable using
%   a symbolic display. 
%
% Examples:
%
%      % Create a scalar variable bounded from below by zero:
%      x = optimvar('x', 'LowerBound', 0);
%      showbounds(x)
%
%           0 <= x
%
%      % Create a 3-by-4 matrix of optimization variables:
%      x = optimvar('x', 3, 4);
%      x(:,1).UpperBound = 10;
%      showbounds(x)
%
%           x(1, 1) <= 10
%           x(2, 1) <= 10
%           x(3, 1) <= 10
%
%      % Create a 3-by-2 matrix of optimization variables with string
%      % indicies and a lower bound of zero:
%      products = {'shirts', 'trousers', 'socks'};
%      locations = {'inside', 'outsourced'};
%      amt = optimvar('amt', products, locations, 'LowerBound', 0);
%      showbounds(amt)
%
%           0 <= amt('shirts', 'inside')
%           0 <= amt('trousers', 'inside')
%           0 <= amt('socks', 'inside')
%           0 <= amt('shirts', 'outsourced')
%           0 <= amt('trousers', 'outsourced')
%           0 <= amt('socks', 'outsourced')
%
%
%   See also OPTIMVAR, SHOWVAR, WRITEBOUNDS, WRITEVAR

 
%   Copyright 2017 The MathWorks, Inc.

