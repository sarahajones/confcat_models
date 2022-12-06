%FINDINDEX Find numeric index in IndexNames
%
%   [NUMIDX1, ... , NUMIDXN] = FINDINDEX(X, STRIDX1, ..., STRIDXN) finds
%   the numeric subscript index equivalent of the specified string indices
%   in the variable X. STRIDXi is either a cellstr, a string vector or a
%   numeric vector. Note that N must match the number of dimensions in X.
%
%   NUMIDX = FINDINDEX(X, STRIDX1, ..., STRIDXN) finds the linear subscript
%   index equivalent of the specified string indices in the variable X.
%   Note that STRIDX1, ..., STRIDXN must have the same size, and NUMIDX
%   will have the same size as STRIDX1, ... STRIDXN. For the array X, if
%   NUMIDX = findindex(X, STRIDX1, ... STRIDXN)) then
%   X(NUMIND(k)) = X(STRIDX1(k), ... STRIDXN(k)) for all k.
%
%   Examples:
%
%   % For repeatability
%   rng(0) 
%
%   % Create and solve an optimization problem
%   p = optimproblem('ObjectiveSense', 'maximize');
%   flow = optimvar('flow', ...
%      {'apples', 'oranges', 'bananas', 'berries'}, {'NYC', 'BOS', 'LAX'}, ...
%      'LowerBound', 0);
%   p.Objective = sum(sum(rand(4, 3).*flow));
%   p.Constraints.NYC = rand(1, 4)*flow(:, 'NYC') <= 10;
%   p.Constraints.BOS = rand(1, 4)*flow(:, 'BOS') <= 12;
%   p.Constraints.LAX = rand(1, 4)*flow(:, 'LAX') <= 35;
%   sol = solve(p);
%
%   % Find the optimal flow of oranges and berries to New York and Los Angeles
%   [idxFruit,idxAirports] = findindex(flow, {'oranges','berries'}, {'NYC', 'LAX'})
%   orangeBerries = sol.flow(idxFruit, idxAirports)
%
%   % List the optimal flow of the following:
%   %  Fruit     Airports
%   %  -----     --------
%   %  Berries   NYC
%   %  Apples    BOS
%   %  Oranges   LAX
%
%   idx = findindex(flow, {'berries', 'apples', 'oranges'}, {'NYC', 'BOS', 'LAX'})
%   optimalFlow = sol.flow(idx)

 
%   Copyright 2017 The MathWorks, Inc.

