function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,x,options,varargin)
%FSOLVE solves systems of nonlinear equations of several variables.
%
%   FSOLVE attempts to solve equations of the form:
%             
%   F(X) = 0    where F and X may be vectors or matrices.   
%
%   FSOLVE implements three different algorithms: trust region dogleg,
%   trust region, and Levenberg-Marquardt. Choose one via the option
%   Algorithm: for instance, to choose trust region, set OPTIONS =
%   optimoptions('fsolve','Algorithm','trust-region'), and then pass
%   OPTIONS to FSOLVE.
%    
%   X = FSOLVE(FUN,X0) starts at the matrix X0 and tries to solve the 
%   equations in FUN.  FUN accepts input X and returns a vector (matrix) of 
%   equation values F evaluated at X. 
%
%   X = FSOLVE(FUN,X0,OPTIONS) solves the equations with the default
%   optimization parameters replaced by values in OPTIONS, an argument
%   created with the OPTIMOPTIONS function.  See OPTIMOPTIONS for details.
%   Use the SpecifyObjectiveGradient option to specify that FUN also
%   returns a second output argument J that is the Jacobian matrix at the
%   point X. If FUN returns a vector F of m components when X has length n,
%   then J is an m-by-n matrix where J(i,j) is the partial derivative of
%   F(i) with respect to x(j). (Note that the Jacobian J is the transpose
%   of the gradient of F.)
%
%   X = FSOLVE(PROBLEM) solves system defined in PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fsolve' in PROBLEM.solver.  Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. 
%
%   [X,FVAL] = FSOLVE(FUN,X0,...) returns the value of the equations FUN 
%   at X. 
%
%   [X,FVAL,EXITFLAG] = FSOLVE(FUN,X0,...) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are listed below. See the documentation
%   for a complete description.
%
%     1  FSOLVE converged to a root.
%     2  Change in X too small.
%     3  Change in residual norm too small.
%     4  Computed search direction too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  Converged to a point that is not a root.
%    -3  Trust region radius too small (Trust-region-dogleg).
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FSOLVE(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, the number of CG iterations (if used) in 
%   OUTPUT.cgiterations, the first-order optimality (if used) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,JACOB] = FSOLVE(FUN,X0,...) returns the 
%   Jacobian of FUN at X.  
%
%   Examples
%     FUN can be specified using @:
%        x = fsolve(@myfun,[2 3 4],optimoptions('fsolve','Display','iter'))
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an anonymous function:
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimoptions('fsolve','Display','off'))
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the system of 
%   nonlinear equations given in the function myfun, which is parameterized 
%   by its second argument c. Here myfun is a MATLAB file function such as
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%           
%   To solve the system of equations for a specific value of c, first 
%   assign the value to c. Then create a one-argument anonymous function 
%   that captures that value of c and calls myfun with two arguments. 
%   Finally, pass this anonymous function to FSOLVE:
%
%       c = -1; % define parameter first
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   See also OPTIMOPTIONS, LSQNONLIN, @.

%   Copyright 1990-2017 The MathWorks, Inc.

% ------------Initialization----------------
defaultopt = struct(...
    'Algorithm','trust-region-dogleg',...
    'DerivativeCheck','off',...
    'Diagnostics','off',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','final',...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward',...
    'FunValCheck','off',...
    'InitDamping', 0.01, ...
    'Jacobian','off',...
    'JacobMult',[],... 
    'JacobPattern','sparse(ones(Jrows,Jcols))',...
    'MaxFunEvals',[],...
    'MaxIter',400,...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))',...
    'OutputFcn',[],...
    'PlotFcns',[],...
    'PrecondBandWidth',Inf,...
    'ScaleProblem','none',...
    'TolFun', 1e-6,...
    'TolFunValue', 1e-6, ...    
    'TolPCG',0.1,...
    'TolX',1e-6,...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel', false );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3
    options = []; 
end

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fsolve:InputArg'));
    end
end

% No options passed. Set options directly to defaultopt after
allDefaultOpts = isempty(options);

% The 'trust-region-reflective' algorithm has been renamed to
% 'trust-region'. Options objects will handle this conversion. However, we
% need to check any structures which may undergo no checking before this
% point.
if isstruct(options) && isfield(options, 'Algorithm') && ...
        ~iscell(options.Algorithm) && ...
        strcmp(options.Algorithm, 'trust-region-reflective')
    options.Algorithm = 'trust-region';
end

% Prepare the options for the solver
[options, optionFeedback] = prepareOptionsForSolver(options, 'fsolve');

% After processing options for optionFeedback, etc., set options to default
% if no options were passed. 
if allDefaultOpts
    % Options are all default
    options = defaultopt;
end

if nargin == 0
    error(message('optim:fsolve:NotEnoughInputs'))
end
% Check for non-double inputs
msg = isoptimargdbl('FSOLVE', {'X0'}, x);
if ~isempty(msg)
    error('optim:fsolve:NonDoubleInput',msg);
end

LB = []; UB = [];
sizes.xShape = size(x);
xstart = x(:);
sizes.nVar = length(xstart);
sizes.mNonlinEq = 0; sizes.mNonlinIneq = 0; % No nonlinear constraints

display = optimget(options,'Display',defaultopt,'fast',allDefaultOpts);
detailedExitMsg = contains(display,'detailed');
switch display
    case {'off','none'}
        verbosity = 0;
    case {'iter','iter-detailed'}
        verbosity = 2;
    case {'final','final-detailed'}
        verbosity = 1;
    case 'testing'
        verbosity = Inf;
    otherwise
        verbosity = 1;
end
diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast',allDefaultOpts),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast',allDefaultOpts),'on');

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast',allDefaultOpts),'on');

algorithm = optimget(options,'Algorithm',defaultopt,'fast',allDefaultOpts);
initDamping = optimget(options,'InitDamping',defaultopt,'fast',allDefaultOpts);
if ~iscell(algorithm)
    initLMparam = initDamping;
else
    initLMparam = algorithm{2}; % Initial Levenberg-Marquardt parameter
    algorithm = algorithm{1};   % Algorithm string
end

switch algorithm
    case 'trust-region-dogleg'
         algorithmflag = 2;
    case 'trust-region'
        algorithmflag = 1;
    case 'levenberg-marquardt'
        algorithmflag = 3;
    otherwise % Invalid choice of Algorithm
        error(message('optim:fsolve:InvalidAlgorithm'))
end

% Process user function
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = lsqfcnchk(FUN,'fsolve',length(varargin),funValCheck,gradflag);
else
    error(message('optim:fsolve:InvalidFUN'))
end

mtxmpy = optimget(options,'JacobMult',defaultopt,'fast',allDefaultOpts);
% Check if name clash
functionNameClashCheck('JacobMult',mtxmpy,'atamult','optim:fsolve:JacobMultNameClash');

% Use internal Jacobian-multiply function if user does not provide JacobMult function 
% or options.Jacobian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungrad') && ~strcmpi(funfcn{1},'fun_then_grad'))
    mtxmpy = @atamult;
end

JAC = [];
x(:) = xstart;
switch funfcn{1}
    case 'fun'
        try
            fuser = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
               getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    case 'fungrad'
        try
            [fuser,JAC] = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
                getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    case 'fun_then_grad'
        try
            fuser = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
                getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
        try
            JAC = feval(funfcn{4},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:JacobianError', ...
                getString(message('optim:fsolve:JacobianError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    otherwise
        error(message('optim:fsolve:UndefinedCalltype'))
end

% Check for non-double data typed values returned by user functions 
if ~isempty( isoptimargdbl('FSOLVE', {'F','J'}, fuser, JAC) )
    error('optim:fsolve:NonDoubleFunVal',getString(message('optimlib:commonMsgs:NonDoubleFunVal','FSOLVE')));
end

f = fuser(:);
sizes.nFun = length(f);

if gradflag
    % check size of JAC
    [Jrows, Jcols] = size(JAC);
    if isempty(options.JacobMult)
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows ~= sizes.nFun || Jcols ~= sizes.nVar
            error(message('optim:fsolve:InvalidJacobian', sizes.nFun, sizes.nVar))
        end
    end
else
    Jrows = sizes.nFun;
    Jcols = sizes.nVar;
end

caller = 'fsolve';

% Choose what algorithm to run: determine algorithmflag and check criteria
if algorithmflag == 1 && sizes.nFun < sizes.nVar
    % trust-region algorithm and not enough equations - switch to
    % levenberg-marquardt algorithm
    warning(message('optim:fsolve:FewerFunsThanVars'))
    algorithmflag = 3;
elseif algorithmflag == 2 && sizes.nFun ~= sizes.nVar
    warning(message('optim:fsolve:NonSquareSystem'));
    algorithmflag = 3;
end

% Set up for diagnostics and derivative check
confcn = {''};
if diagnostics
    % Do diagnostics on information so far
    constflag = false; gradconstflag = false; hessflag = false; 
    non_eq = 0;non_ineq = 0;lin_eq = 0;lin_ineq = 0;
   
    % Set OUTPUT.algorithm for diagnostics
     switch algorithmflag
         case 1
             OUTPUT.algorithm = 'trust-region';
         case 2
             OUTPUT.algorithm = 'trust-region-dogleg';
         case 3
             OUTPUT.algorithm = 'levenberg-marquardt';
    end
    diagnose('fsolve',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
        xstart,non_eq,non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn);
end

% Read options for finitedifferences
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast',allDefaultOpts);
options.GradObj = optimget(options,'Jacobian',defaultopt,'fast',allDefaultOpts);
options.GradConstr = 'off';
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast',allDefaultOpts),'on');
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx;
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast',allDefaultOpts);
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast',allDefaultOpts);
options = validateFinDiffRelStep(sizes.nVar,options,defaultopt);
options.UseParallel = optimget(options,'UseParallel',defaultopt,'fast',allDefaultOpts);

% Create structure of flags for finitedifferences
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward'); % Check for forward fin-diff
finDiffFlags.scaleObjConstr = false; % No scaling
finDiffFlags.chkComplexObj = false;  % Don't check whether objective function values are complex
finDiffFlags.isGrad = false;         % Compute Jacobian, not gradient
finDiffFlags.hasLBs = false(sizes.nVar,1);  % No bounds
finDiffFlags.hasUBs = false(sizes.nVar,1);
% finDiffFlags.chkFunEval will be set below according to 
% where finDiffFlags.chkFunEval will be used.

% Check derivatives
if DerivativeCheck && gradflag          % user wants to check derivatives
    lb = -Inf(sizes.nVar,1); ub = Inf(sizes.nVar,1);
    finDiffFlags.chkFunEval = false;    % don't check function values during finite differences
    validateFirstDerivatives(funfcn,confcn,x, ...
        lb,ub,options,finDiffFlags,sizes,varargin{:});
end

% Flag to determine whether to look up the exit msg.
makeExitMsg = logical(verbosity) || nargout > 3;

% Execute algorithm
if algorithmflag == 1   % trust-region
    if ~gradflag
        Jstr = optimget(options,'JacobPattern',defaultopt,'fast',allDefaultOpts);
        if ischar(Jstr)
            % options.JacobPattern is the default: 'sparse(ones(jrows,jcols))'
            % Bypass column re-ordering, graph coloring, etc. since we will
            % perform dense finite-differencing anyway.
            Jstr = [];
            finDiffFlags.sparseFinDiff = false;
        else
            checkoptionsize('JacobPattern', size(Jstr), Jcols, Jrows);
            finDiffFlags.sparseFinDiff = true;
        end
    else
        Jstr = [];
    end
    computeLambda = 0;
    % Set MaxFunEvals appropriately for trust-region
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    % Don't check function values during finite differences
    finDiffFlags.chkFunEval = false;     
    [x,FVAL,~,JACOB,EXITFLAG,OUTPUT,msgData]=...
        snls(funfcn,x,LB,UB,verbosity,options,defaultopt,f,JAC,caller,Jstr,...
        computeLambda,mtxmpy,detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
    
    % Correct the algorithm name stored in snls
    OUTPUT.algorithm = 'trust-region';
elseif algorithmflag == 2   % trust-region dogleg    
    % Set MaxFunEvals appropriately for trust-region-dogleg
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    % Check function values during finite differences
    finDiffFlags.chkFunEval = true;    
    
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msgData]=...
        trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,f,JAC,...
        detailedExitMsg,optionFeedback,finDiffFlags,sizes,varargin{:});
elseif algorithmflag == 3   % Levenberg-Marquardt
    % Set MaxFunEvals appropriately for LM
    defaultopt.MaxFunEvals = '200*numberOfVariables';
    % Check function values during finite differences
    finDiffFlags.chkFunEval = true;    
    
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msgData] = ...
        levenbergMarquardt(funfcn,x,verbosity,options,defaultopt,f,JAC,caller, ...
        initLMparam,detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
end

Resnorm = FVAL'*FVAL;  % assumes FVAL still a vector
sqrtTolFunValue = sqrt(optimget(options,'TolFunValue',defaultopt,'fast',allDefaultOpts));
if EXITFLAG > 0 % if we think we converged:
    % Call createExitMsg with appended additional information on the closeness
    % to a root.
    if Resnorm > sqrtTolFunValue
        msgData = internalFlagForExitMessage(algorithmflag == 2,msgData,EXITFLAG);
        EXITFLAG = -2;
    end   
    if makeExitMsg
        OUTPUT.message = createExitMsg(msgData{:},Resnorm,optionFeedback.TolFunValue,sqrtTolFunValue);
    end
else
    if makeExitMsg
        OUTPUT.message = createExitMsg(msgData{:});
    end
end

% Reset FVAL to shape of the user-function output, fuser
FVAL = reshape(FVAL,size(fuser));
end

%--------------------------------------------------------------------------
function msgData = internalFlagForExitMessage(isDogleg,msgData,EXITFLAG)
% internalFlagForExitMessage Create internal flag for exit message.
% Change internal exitflag to unique identifier -21, -22, -23 or -27

if isDogleg && msgData{2} == 26
    % Algorithm dogleg, stopped because of undefined Jacobian
    msgData{2} = -27;
else
    % All algorithms: change internal exitflag by negating the
    % exitflag and adding to -20.
    msgData{2} = -20 - EXITFLAG;
end
end
