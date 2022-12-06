function W = hmfleq1(Hinfo,Y,V)
%HMFLEQ1 Hessian-matrix product function for BROWNVV objective.
% Documentation example.
% W = hmfbx4(Hinfo,Y,V) computes W = (Hinfo-V*V')*Y
% where Hinfo is a sparse matrix computed by BROWNVV 
% and V is a 2 column matrix.

%   Copyright 1984-2008 The MathWorks, Inc.
  
W = Hinfo*Y - V*(V'*Y);
  



