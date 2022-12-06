%SHOWEXPR Display an OptimizationExpression
%
%  SHOWEXPR(EXPR) displays a simplified mathematical form of the
%  OptimizationExpression EXPR.
%
%  For example,
%  amt = optimvar('amt',{'shirts','trousers','socks'}, ...
%                       {'inside','outsourced'});
%
%  unitcost = [0.6, 0.8;
%              0.8, 0.9;
%              0.3, 0.4];
%
%  totalCostPerItem = sum(unitcost.*amt, 2);
%
%  showexpr(totalCostPerItem)
%
%  (1, 1)
%
%    0.6*amt('shirts', 'inside') + 0.8*amt('shirts', 'outsourced')
%
%  (2, 1)
%
%    0.8*amt('trousers', 'inside') + 0.9*amt('trousers', 'outsourced')
%
%  (3, 1)
%
%    0.3*amt('socks', 'inside') + 0.4*amt('socks', 'outsourced')
%
%  See also OPTIMEXPR, WRITEEXPR

 
%   Copyright 2017 The MathWorks, Inc.

