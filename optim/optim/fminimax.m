function [x,FVAL,MAXFVAL,EXITFLAG,OUTPUT,LAMBDA] = fminimax(FUN,x,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FMINIMAX finds a minimax solution of a function of several variables. 
%   FMINIMAX attempts to solve the following problem:
%   min (max {FUN(X} )  where FUN and X can be vectors or matrices.
%    X 
% 
%   X = FMINIMAX(FUN,X0) starts at X0 and finds a minimax solution X to 
%   the functions in FUN. FUN accepts input X and returns a vector
%   (matrix) of function values F evaluated at X. X0 may be a scalar,
%   vector, or matrix. 
%
%   X = FMINIMAX(FUN,X0,A,B) solves the minimax problem subject to the
%   linear inequalities A*X <= B.
%
%   X = FMINIMAX(FUN,X0,A,B,Aeq,Beq) solves the minimax problem
%   subject to the linear equalities Aeq*X = Beq as well.  (Set A = [] and 
%   B = [] if no inequalities exist.)
%
%   X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower 
%   and upper bounds on the design variables, X, so that the solution is 
%   in the range LB <= X <= UB. You may use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%   set UB(i) = Inf if X(i) is unbounded above.
%   
%   X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) subjects the 
%   goal attainment problem to the constraints defined in NONLCON (usually 
%   a MATLAB file: NONLCON.m). The function NONLCON should return the vectors
%   C and Ceq, representing the nonlinear inequalities and equalities 
%   respectively, when called with feval: [C, Ceq] = feval(NONLCON,X). 
%   FMINIMAX optimizes such that C(X) <= 0 and Ceq(X) = 0.
%
%   X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
%   the default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details. Use the SpecifyObjectiveGradient option to specify that FUN
%   may be called with two output arguments where the second, G, is the
%   partial derivatives of the function df/dX, at the point X: [F,G] =
%   feval(FUN,X). Use the SpecifyConstraintGradient option to specify that
%   NONLCON may be called with four output arguments: [C,Ceq,GC,GCeq] =
%   feval(NONLCON,X) where GC is the partial derivatives of  the constraint
%   vector of inequalities C an GCeq is the partial derivatives  of the
%   constraint vector of equalities Ceq. Use OPTIONS = [] as a place
%   holder if no options are set.
%
%   X = FMINIMAX(PROBLEM) finds a minimax solution for PROBLEM. PROBLEM is 
%   a structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the linear inequality constraints in PROBLEM.Aineq
%   and PROBLEM.bineq, the linear equality constraints in PROBLEM.Aeq and
%   PROBLEM.beq, the lower bounds in PROBLEM.lb, the upper bounds in 
%   PROBLEM.ub, the nonlinear constraint function in PROBLEM.nonlcon, the
%   options structure in PROBLEM.options, and solver name 'fminimax' in
%   PROBLEM.solver. Use this syntax to solve at the command line a problem 
%   exported from OPTIMTOOL. 
%
%   [X,FVAL] = FMINIMAX(FUN,X0,...) returns the value of the objective 
%   functions at the solution X: FVAL = feval(FUN,X).
%
%   [X,FVAL,MAXFVAL] = FMINIMAX(FUN,X0,...) returns 
%   MAXFVAL = max { FUN(X) } at the solution X.
%
%   [X,FVAL,MAXFVAL,EXITFLAG] = FMINIMAX(FUN,X0,...) returns an EXITFLAG
%   that describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are listed below. See the documentation
%   for a complete description.
%
%     1  FMINIMAX converged to a solution.
%     4  Computed search direction too small.
%     5  Predicted change in max objective function too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  No feasible point found.
%   
%   [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT] = FMINIMAX(FUN,X0,...) returns a 
%   structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the norm of the final step in OUTPUT.stepsize, the 
%   final line search steplength in OUTPUT.lssteplength, the algorithm used
%   in OUTPUT.algorithm, the first-order optimality in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message. 
%
%   [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINIMAX(FUN,X0,...) returns 
%   the Lagrange multipliers at the solution X: LAMBDA.lower for LB, 
%   LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear inequalities, 
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear 
%   equalities.
%
%   Examples
%     FUN can be specified using @:
%        x = fminimax(@myfun,[2 3 4])
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = cos(x);
%
%   FUN can also be an anonymous function:
%
%       x = fminimax(@(x) sin(3*x),[2 5])
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve a minimax 
%   problem where the objectives given in the function myfun are 
%   parameterized by its second argument c. Here myfun is a MATLAB file 
%   function such as
%
%       function F = myfun(x,c)
%       F = [x(1)^2 + c*x(2)^2;
%            x(2) - x(1)];
%
%   To optimize for a specific value of c, first assign the value to c. 
%   Then create a one-argument anonymous function that captures that value 
%   of c and calls myfun with two arguments. Finally pass this anonymous 
%   function to FMINIMAX:
%
%       c = 2; % define parameter first
%       x = fminimax(@(x) myfun(x,c),[1;1])
%
%   See also OPTIMOPTIONS, @, INLINE, FGOALATTAIN, LSQNONLIN.

%   Copyright 1990-2017 The MathWorks, Inc.

defaultopt = struct( ...
    'Diagnostics','off', ...
    'DiffMaxChange',Inf, ...
    'DiffMinChange',0, ...
    'Display','final', ...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward', ...
    'FunValCheck','off', ...
    'GradConstr','off', ...
    'GradObj','off', ...
    'Hessian','off', ...    % Not used
    'LargeScale','off', ... % Not used    
    'MaxFunEvals','100*numberOfVariables', ...
    'MaxIter',400, ...    
    'MaxSQPIter','10*max(numberOfVariables,numberOfInequalities+numberOfBounds)', ...    
    'MeritFunction','multiobj', ...    
    'MinAbsMax',0, ...
    'OutputFcn',[], ...
    'PlotFcns',[], ...    
    'RelLineSrchBnd',[], ...
    'RelLineSrchBndDuration',1, ...    
    'TolCon',1e-6, ...
    'TolConSQP',1e-6, ...
    'TolFun',1e-6, ...
    'TolFunValue',1e-6, ...
    'TolX',1e-6, ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel',false ...
    );
% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 10
    options = [];
    if nargin < 9
        NONLCON = [];
        if nargin < 8
            UB = [];
            if nargin < 7
                LB = [];
                if nargin < 6
                    Beq = [];
                    if nargin < 5
                        Aeq = [];
                        if nargin < 4
                            B = [];
                            if nargin < 3
                                A = [];
                            end
                        end
                    end
                end
            end
        end
    end
end

algAS = 'active-set';

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,A,B,Aeq,Beq,LB,UB,NONLCON,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fminimax:InputArg'));
    end
end

if nargin == 0 
    error(message('optim:fminimax:NotEnoughInputs'))
end

% No options passed. Set options directly to defaultopt after
allDefaultOpts = isempty(options);

% Prepare the options for the solver
[options, optionFeedback] = prepareOptionsForSolver(options, 'fminimax');

% Check for non-double inputs
msg = isoptimargdbl('FMINIMAX', {'X0','A','B','Aeq','Beq','LB','UB'}, ...
                                  x,   A,  B,  Aeq,  Beq,  LB,  UB);
if ~isempty(msg)
    error('optim:fminimax:NonDoubleInput',msg);
end

% After processing options for optionFeedback, etc., set options to default
% if no options were passed. 
if allDefaultOpts
    % Options are all default
    options = defaultopt;
end

initVals.xOrigShape = x;
xnew = [x(:); 0];

sizes.xShape = size(x);
sizes.nVar = length(x(:));
numberOfVariablesplus1 = length(xnew);

diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast',allDefaultOpts),'on');

display = optimget(options,'Display',defaultopt,'fast',allDefaultOpts);
flags.detailedExitMsg = contains(display,'detailed');
switch display
    case {'off','none'}
        verbosity = 0;
    case {'notify','notify-detailed'}
        verbosity = 1;
    case {'final','final-detailed'}
        verbosity = 2;
    case {'iter','iter-detailed'}
        verbosity = 3;
    otherwise
        verbosity = 2;
end

% Set to column vectors
B = B(:);
Beq = Beq(:);

[xnew(1:sizes.nVar),l,u,msg] = checkbounds(xnew(1:sizes.nVar),LB,UB,sizes.nVar);
if ~isempty(msg)
    EXITFLAG = -2;
    [FVAL,MAXFVAL,LAMBDA] = deal([]);
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.stepsize = [];
    OUTPUT.lssteplength = [];
    OUTPUT.algorithm = algAS;
    OUTPUT.firstorderopt = [];
    OUTPUT.constrviolation =[];
    OUTPUT.message = msg;
    x(:) = xnew(1:sizes.nVar);
    if verbosity > 0
        disp(msg)
    end
    return
end

neqgoals = optimget(options, 'MinAbsMax',defaultopt,'fast',allDefaultOpts);

% flags.meritFunction is 1 unless changed by user to fmincon merit function;
% formerly options(7)
% 0 uses the fmincon single-objective merit and Hess; 1 is the default
flags.meritFunction = strcmp(optimget(options,'MeritFunction',defaultopt,'fast',allDefaultOpts),'multiobj');
lenVarIn = length(varargin);
% goalcon and goalfun also take:
% neqgoals,funfcn,gradfcn,WEIGHT,GOAL,x,errCheck
goalargs = 7; 

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast',allDefaultOpts),'on');
% Gather options needed for finitedifferences
% Write checked DiffMaxChange, DiffMinChage, FinDiffType, FinDiffRelStep,
% GradObj and GradConstr options back into struct for later use
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast',allDefaultOpts);
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast',allDefaultOpts);
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast',allDefaultOpts);
if options.DiffMinChange >= options.DiffMaxChange
    error(message('optim:fgoalattain:DiffChangesInconsistent', sprintf( '%0.5g', options.DiffMinChange ), sprintf( '%0.5g', options.DiffMaxChange )))
end
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx(:);
options = validateFinDiffRelStep(sizes.nVar,options,defaultopt);

options.GradObj = optimget(options,'GradObj',defaultopt,'fast',allDefaultOpts);
options.GradConstr = optimget(options,'GradConstr',defaultopt,'fast',allDefaultOpts);

flags.grad = strcmp(options.GradObj,'on');
flags.gradconst = strcmp(options.GradConstr,'on');
if strcmpi(optimget(options,'Hessian',defaultopt,'fast',allDefaultOpts),'on')
    warning(message('optim:fminimax:UserHessNotUsed'))
end
flags.hess = false; % Algorithm uses BFGS Hessian approximation
userconstflag = ~isempty(NONLCON);

% If nonlinear constraints exist, need either both function and constraint
% gradients, or none
if userconstflag
    flags.gradconst = flags.grad && flags.gradconst;
else % No user nonlinear constraints
    flags.gradconst = flags.grad;
end
flags.grad = true; % Always can compute gradient of goalfun since based on x

% Update options GradObj and GradConstr to reflect the update for the
% constraint function
if ~flags.gradconst
    options.GradObj = 'off';
    options.GradConstr = 'off';
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    % Pass flags.gradconst as the flag which tells whether or not to
    % evaluate gradients from the user function. flags.grad is meant for
    % goalfun and is always set to true for this problem.
   funfcn = optimfcnchk(FUN,'goalcon',length(varargin),funValCheck,flags.gradconst,flags.hess);
else
   error(message('optim:fminimax:invalidFUN'))
end
% We can always compute gradient since based only on xnew.
% Pass in false for funValCheck argument as goalfun is not a user function.
ffun = optimfcnchk(@goalfun,'fminimax',lenVarIn+goalargs,false,flags.grad);

if userconstflag % NONLCON is non-empty, goalcon is the caller to NONLCON
   confcn = ...
      optimfcnchk(NONLCON,'goalcon',length(varargin),funValCheck,flags.gradconst,false,true);
else
   confcn{1} = '';
end
% Pass in false for funValCheck argument as goalfun is not a user function
cfun = optimfcnchk(@goalcon,'fminimax',lenVarIn+goalargs,false,flags.gradconst,false,true); 

lenvlb = length(l);
lenvub = length(u);

i = 1:lenvlb;
lindex = xnew(i) < l(i);
if any(lindex)
   xnew(lindex) = l(lindex) + 1e-4; 
end
i = 1:lenvub;
uindex = xnew(i) > u(i);
if any(uindex)
   xnew(uindex) = u(uindex);
end
x(:) = xnew(1:end-1);

% Evaluate user function to get number of function values at x.
user_f = feval(funfcn{3},x,varargin{:});
user_f = user_f(:);
sizes.nFun = length(user_f);

% Check if neqgoals (MinAbsMax) is less or equal to the length of user function                           
if neqgoals > sizes.nFun
    warning(message('optim:fminimax:InconsistentNumEqGoal'))
    % The number of F(x) to minimize the worst case absolute values can be
    % at most equal to the length of user objective function.
    neqgoals = sizes.nFun;
end

WEIGHT = ones(sizes.nFun,1);
GOAL = zeros(sizes.nFun,1);

initVals.g = zeros(numberOfVariablesplus1,1);
initVals.H = [];
errCheck = true; % Perform error checking on initial function evaluations

extravarargin = [{neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck}, varargin]; 
% Evaluate goal function
switch ffun{1}
    case 'fun'
        initVals.f = feval(ffun{3},xnew,extravarargin{:});
    case 'fungrad'
        [initVals.f,initVals.g] = feval(ffun{3},xnew,extravarargin{:});
    otherwise
        error(message('optim:fminimax:UndefinedCalltype'))
end

% Evaluate goal constraints
switch cfun{1}
    case 'fun'
        [ctmp,ceqtmp] = feval(cfun{3},xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
        initVals.gnc = zeros(numberOfVariablesplus1,length(initVals.ncineq));
        initVals.gnceq = zeros(numberOfVariablesplus1,length(initVals.nceq));
    case 'fungrad'
        [ctmp,ceqtmp,initVals.gnc,initVals.gnceq] = feval(cfun{3},xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
    otherwise
        error(message('optim:fminimax:UndefinedCalltype'))
end

% Make sure empty constraint and their derivatives have correct sizes (not 0-by-0):
if isempty(initVals.ncineq)
    initVals.ncineq = reshape(initVals.ncineq,0,1);
end
if isempty(initVals.nceq)
    initVals.nceq = reshape(initVals.nceq,0,1);
end
if isempty(Aeq)
    Aeq = reshape(Aeq,0,sizes.nVar);
    Beq = reshape(Beq,0,1);
end
if isempty(A)
    A = reshape(A,0,sizes.nVar);
    B = reshape(B,0,1);    
end

sizes.mNonlinEq = length(initVals.nceq);
sizes.mNonlinIneq = length(initVals.ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);

if Aeqcol ~= sizes.nVar
   error(message('optim:fminimax:InvalidSizeOfAeq', sizes.nVar))
end
if Acol ~= sizes.nVar
   error(message('optim:fminimax:InvalidSizeOfA', sizes.nVar))
end

just_user_constraints = sizes.mNonlinIneq - sizes.nFun - neqgoals;
OUTPUT.algorithm = algAS;

if diagnostics
    % Do diagnostics on information so far
    diagnose('fminimax',OUTPUT,flags.gradconst,flags.hess,userconstflag,flags.gradconst,...
        xnew(1:end-1),sizes.mNonlinEq,just_user_constraints,lin_eq,lin_ineq,l,u,funfcn,confcn);
end

% Add extra column to account for extra xnew component
A = [A,zeros(lin_ineq,1)];
Aeq = [Aeq,zeros(lin_eq,1)];

% Only need to perform error checking on initial function evaluations
errCheck = false;

% Convert function handles to anonymous functions with additional arguments
% in its workspace. Even though ffun and cfun are internal functions, put fevals
% here for consistency.
ffun{3} = @(y,varargin) feval(ffun{3},y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});
cfun{3} = @(y,varargin) feval(cfun{3},y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});

% Problem related data is passed to nlconst in problemInfo structure
problemInfo.nHardConstraints = neqgoals;
problemInfo.weight = WEIGHT;
problemInfo.goal = GOAL;

% Create default structure of flags for finitedifferences:
% This structure will (temporarily) ignore some of the features that are
% algorithm-specific (e.g. scaling and fault-tolerance) and can be turned
% on later for the main algorithm.
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward');
finDiffFlags.scaleObjConstr = false; % No scaling for now
finDiffFlags.chkFunEval = false;     % No fault-tolerance yet
finDiffFlags.chkComplexObj = false;  % No need to check for complex values
finDiffFlags.isGrad = false;         % Multi-objective
finDiffFlags.hasLBs = false(sizes.nVar,1);
finDiffFlags.hasUBs = false(sizes.nVar,1);
if ~isempty(l)
    finDiffFlags.hasLBs = isfinite(l);   % Finite lower bounds
end
if ~isempty(u)
    finDiffFlags.hasUBs = isfinite(u);   % Finite upper bounds
end

% Adjust variables used by finite-differencing for problem re-formulation:
% i.e. adjust nVar-length vectors for auxiliary variable
options.TypicalX = [typicalx(:); 1];
if finDiffFlags.fwdFinDiff
    options.FinDiffRelStep = [options.FinDiffRelStep; sqrt(eps)];
else
    options.FinDiffRelStep = [options.FinDiffRelStep; eps^(1/3)];
end
l = [l;-Inf];
u = [u; Inf];
finDiffFlags.hasLBs = [finDiffFlags.hasLBs; false];
finDiffFlags.hasUBs = [finDiffFlags.hasUBs; false];
finDiffFlags.isGrad = true;         % New formulation has single objective

% For parallel finite difference (if needed) we need to send the function
% handles now to the workers. This avoids sending the function handles in
% every iteration of the solver. The output from 'setOptimFcnHandleOnWorkers' 
% is a onCleanup object that will perform cleanup task on the workers.
UseParallel = optimget(options,'UseParallel',defaultopt,'fast',allDefaultOpts);
cleanupObj = setOptimFcnHandleOnWorkers(UseParallel,ffun,cfun); %#ok<NASGU>

% Flag to determine whether to look up the exit msg.
flags.makeExitMsg = logical(verbosity) || nargout > 4;

[xnew,~,LAMBDA,EXITFLAG,OUTPUT]=...
   nlconst(ffun,xnew,l,u,full(A),B,full(Aeq),Beq,cfun,options,defaultopt, ...
   finDiffFlags,verbosity,flags,initVals,problemInfo,optionFeedback,varargin{:});

if ~isempty(LAMBDA)
    just_user_constraints = length(LAMBDA.ineqnonlin) - sizes.nFun - neqgoals;
    LAMBDA.ineqnonlin = LAMBDA.ineqnonlin(1:just_user_constraints);
    LAMBDA.lower = LAMBDA.lower(1:sizes.nVar);
    LAMBDA.upper = LAMBDA.upper(1:sizes.nVar);
end

% Evaluate user objective functions, find maxfval
x(:) = xnew(1:end-1);
FVAL = feval(funfcn{3},x,varargin{:});
MAXFVAL = max(FVAL);

% Force a cleanup of the handle object. Sometimes, MATLAB may
% delay the cleanup but we want to be sure it is cleaned up.
clear cleanupObj
