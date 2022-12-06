%WRITEBOUNDS Write bounds on optimization variables to file
%
%   WRITEBOUNDS(X) writes the finite bounds for OptimizationVariable array
%   X into a text file. The file name is the NAME property of X, appended
%   with '_bounds.txt'. writebounds overwrites any existing file.
%
%   WRITEBOUNDS(X, FILENAME) writes the finite bounds for the variables in
%   X to FILENAME using a symbolic display.
%
% Examples:
%
%      % Create a scalar variable bounded from below by zero:
%      x = optimvar('x', 'LowerBound', 0);
%      writebounds(x, 'text.txt')
%
%      Overwrites (or creates if nonexistent) the local text file
%      'text.txt' with the following text:
%
%           0 <= x
%
%
%      % Create a 3-by-4 matrix of optimization variables:
%      x = optimvar('x', 3, 4);
%      x(:,1).UpperBound = 10;
%      writebounds(x, 'text.txt')
%
%      Overwrites (or creates if nonexistent) the local text file
%      'text.txt' with the following text:
%
%           x(1,1) <= 10
%           x(2,1) <= 10
%           x(3,1) <= 10
%
%
%      % Create a 3-by-2 matrix of optimization variables with string
%      % indicies and a lower bound of zero:
%      products = {'shirts', 'trousers', 'socks'};
%      locations = {'inside', 'outsourced'};
%      amt = optimvar('amt', products, locations, 'LowerBound', 0);
%      writebounds(amt, 'text.txt')
%
%      Overwrites (or creates if nonexistent) the local text file
%      'text.txt' with the following text:
%
%           0 <= amt('shirts','inside')
%           0 <= amt('trousers ','inside')
%           0 <= amt('socks ','inside')
%           0 <= amt('shirts',' outsourced ')
%           0 <= amt('trousers ',' outsourced ')
%           0 <= amt('socks ',' outsourced ')
%
%
%   See also OPTIMVAR, SHOWBOUNDS, SHOWVAR, WRITEVAR

 
%   Copyright 2017-2018 The MathWorks, Inc.

