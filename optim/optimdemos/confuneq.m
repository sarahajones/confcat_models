function [c, ceq] = confuneq(x)
%CONFUNEQ Nonlinear inequality and equality constraints.
% Documentation example.

%   Copyright 1990-2008 The MathWorks, Inc.

c = -x(1)*x(2) - 10;
% Nonlinear equality constraint:
ceq = x(1)^2 + x(2) - 1;
