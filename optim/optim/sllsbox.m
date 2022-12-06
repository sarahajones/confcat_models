function [x,residnorm,residual,fnrm,it,npcg,exitflag,LAMBDA,msg] = ...
    sllsbox(A,b,lb,ub,xstart,params,options,defaultopt,mtxmpy,computeLambda,varargin)
%

%SLLSBOX Linear least-squares with bounds
%
% x=slls(A,b,lb,ub) returns the solution to the
% box-constrained linear least-squares problem,
%
%        min { 0.5*||Ax - b||^2 : lb <= x <= ub}.
%
% where A is a matrix with more rows than cols. (may be virtual)
%

%   Copyright 1990-2018 The MathWorks, Inc.

if nargin < 4
  error(message('optim:sllsbox:NotEnoughInputs'))
end
if nargin < 5
    xstart = []; 
end
if nargin < 6 
    params.verb = 1;
    params.emptyInitialPoint = true;
end 
if nargin < 7
    options = []; 
end

pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tolfun = optimget(options,'TolFun',defaultopt,'fast') ;
tolfunvalue = optimget(options,'TolFunValue',defaultopt,'fast');
itb = optimget(options,'MaxIter',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast');

% Add Preconditioner to defaultopt. This is not a documented option for
% lsqlin, and is not added to defaultopt at the user-facing function level.
defaultopt.Preconditioner = @aprecon;
pcmtx = optimget(options,'Preconditioner',defaultopt,'fast'); 

c = feval(mtxmpy,A,-b,-1,varargin{:}); 
n = length(c); 
it = 1; 
cvec = c; 
% In case the defaults were gathered from calling: optimset('lsqlin'):
numberOfVariables = n;
kmax = optimget(options,'MaxPCGIter',defaultopt,'fast') ;
typx = optimget(options,'TypicalX',defaultopt,'fast') ;
if ischar(kmax)
   if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
      kmax = max(1,floor(numberOfVariables/2));
   else
      error(message('optim:sllsbox:InvalidMaxPCGIter'))
   end
end
if ischar(typx)
   if isequal(lower(typx),'ones(numberofvariables,1)')
      typx = ones(numberOfVariables,1);
   else
      error(message('optim:sllsbox:InvalidTypicalX'))
   end
end
checkoptionsize('TypicalX', size(typx), numberOfVariables);

if n == 0
   error(message('optim:sllsbox:InvalidN'))
end
if isempty(lb), lb = -inf*ones(n,1); end
if isempty(ub), ub = inf*ones(n,1); end
arg = (ub >= 1e10); arg2 = (lb <= -1e10);
ub(arg) = inf;
lb(arg2) = -inf;
if any(ub == lb) 
   error(message('optim:sllsbox:EqualLowerUpperBnd'))
elseif min(ub-lb) <= 0
   error(message('optim:sllsbox:InconsistentBnds'))
end
lvec = lb; uvec = ub;

% Initial point
xinitOutOfBounds_idx = xstart < lb | xstart > ub;
if any(xinitOutOfBounds_idx)
    if params.emptyInitialPoint
        % When no x0 provided, x0 is set to the zero vector in lsqlin.
        % If zero start point violates at least one bound, set x0
        % to box-centered point. (This is done this way to maintain
        % backwards compatibility.)
        xstart = startx(ub,lb);
    else
        % If user-provided x0 has components not within bounds,
        % set those components to a box-centered point
        xstart = startx(ub,lb,xstart,xinitOutOfBounds_idx);
    end
end

if isempty(typx), typx = ones(n,1); end
if isempty(params.verb), params.verb = 1; end

%   SHIFT AND SCALE
[xstart,lb,ub,~,DS,c] = shiftsc(xstart,lb,ub,typx,'sllsbox',mtxmpy,cvec,A,varargin{:});

%   MORE INITIALIZATIONS
tol2 = sqrt(tolfunvalue);
dellow = 1.; delup = 10^3; npcg = 0; digits = inf; 
done = false;
del = 10*eps; posdef = 0;x = xstart; y = x;sigma = ones(n,1);
oval = inf; [val,g] = fquad(x,c,A,'sllsbox',mtxmpy,DS,varargin{:});
prev_diff = 0;
printMsg = params.verb > 0;

%   MAIN LOOP: GENERATE FEAS. SEQ.  x(it) S.T. q(x(it)) IS DECREASING.
while ~done
   %     Update and display
   [v,dv] = definev(g,x,lb,ub);
   fnrm = norm(v.*g,inf); 
   %
   %     TEST FOR CONVERGENCE
   diff = abs(oval-val);
   if it > 1
       digits = prev_diff/max(diff,eps); 
   end
   
   prev_diff = diff; oval = val; 
   if ((fnrm < tolfun) && (posdef == 1)) 
      exitflag = 1; done = true;
      msg = qpTrustRegionMsg('optim:sllsbox:ExitSolved','qp_bound_optimal',printMsg,'');
   elseif diff < tolfunvalue*(1+abs(oval))
      exitflag = 3; done = true;
      msg = qpTrustRegionMsg('optim:sllsbox:ExitDeltaF','qp_functiontolerance',printMsg,'lsqlin');  
   elseif (  (diff < tol2*(1+abs(oval))) && (digits < 3.5))
      exitflag = 3; done = true;
      msg = qpTrustRegionMsg('optim:sllsbox:ExitDeltaFSlowChg','qp_functiontolerance_sqrt',printMsg,'lsqlin');
   end

   %
   if ~done
      %       DETERMINE THE SEARCH DIRECTION
      dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd).*sigma));
      grad = D*g; normg = norm(grad); 
      delta = max(dellow,norm(v)); delta = min(delta,delup);
      [s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,...
         pcmtx,pcflags,pcgtol,A,1,kmax,varargin{:});
      npcg = npcg + pcgit;
      
      %       DO A REFLECTIVE (BISECTION) LINE SEARCH. UPDATE x,y,sigma.
      strg= s'*(sigma.*g); ox = x;  osig = sigma; ostrg = strg;
      if strg >= 0 
         exitflag = -4; done = true;
         msg = qpTrustRegionMsg('optim:sllsbox:ExitIllCond','qp_ill_conditioned',printMsg,'lsqlin');
         if params.verb > 0
            fprintf('%s\n', msg);
         end
      else
         [x,sigma,alpha] = biqpbox(s,c,ostrg,ox,y,osig,lb,ub,oval,posdef,...
            normg,DS,mtxmpy,A,1,varargin{:});
         y = y + alpha*s; 
         
         %          PERTURB x AND y ?
         [~,x,y] = perturbTrustRegionReflective(x,lb,ub,del,y,sigma);
         
         %          EVALUATE NEW FUNCTION VALUE, GRADIENT. 
         it = it + 1; 
         [val,g] = fquad(x,c,A,'sllsbox',mtxmpy,DS,varargin{:}); 
      end
   end
   if it >= itb 
      it = it - 1; 
      exitflag = 0; done = true;
      msg = qpTrustRegionMsg('optim:sllsbox:ExitMaxIter','solver_stopped_premature',printMsg,'lsqlin');
   end
end

%   RESCALE, UNSHIFT, AND EXIT.
x = unshsca(x,lvec,uvec,DS);

residual = feval(mtxmpy,A,x,1,varargin{:}) - b; % A*x-b
residnorm = sum(residual.*residual);

g = feval(mtxmpy,A,residual,-1,varargin{:}); % A'*(A*x-b) so residual
g = full(g);

if computeLambda
    LAMBDA = formBoxLambda(g,lvec,uvec);
    LAMBDA.ineqlin = []; 
    LAMBDA.eqlin = [];
else
    LAMBDA = [];
end
