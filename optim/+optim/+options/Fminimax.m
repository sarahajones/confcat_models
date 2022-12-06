classdef (Sealed) Fminimax < optim.options.SingleAlgorithm
%

%Fminimax Options for FMINIMAX
%
%   The OPTIM.OPTIONS.FMINIMAX class allows the user to create a set of
%   options for the FMINIMAX solver. For a list of options that can be set,
%   see the documentation for FMINIMAX.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX creates a set of options for FMINIMAX
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX(PARAM, VAL, ...) creates a set of options
%   for FMINIMAX with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.SINGLEALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

%   Copyright 2012-2017 The MathWorks, Inc.

    properties (Dependent)

%ABSOLUTEMAXOBJECTIVECOUNT Number of elements of Fi(x) to minimize the
%                          maximum absolute value of Fi
% 
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.                        
        AbsoluteMaxObjectiveCount

%CONSTRAINTTOLERANCE Tolerance on the constraint violation
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        ConstraintTolerance
        
%DISPLAY Level of display
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.        
        Display
        
%FINITEDIFFERENCESTEPSIZE Scalar or vector step size factor
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        FiniteDifferenceStepSize
        
%FINITEDIFFERENCETYPE Finite difference type
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        FiniteDifferenceType

%FUNCTIONTOLERANCE Termination tolerance on the change in function value
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        FunctionTolerance
                        
%MAXITERATIONS Maximum number of iterations allowed
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        MaxIterations
        
%MAXFUNCTIONEVALUATIONS Maximum number of function evaluations allowed
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        MaxFunctionEvaluations

%OPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
%                    measure
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        OptimalityTolerance               
        
%OUTPUTFCN Callbacks that are called at each iteration
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        OutputFcn
        
%PLOTFCN Plots various measures of progress while the algorithm executes
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        PlotFcn

%SPECIFYCONSTRAINTGRADIENT Gradient for nonlinear constraint functions
%                          defined by the user
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        SpecifyConstraintGradient
        
%SPECIFYOBJECTIVEGRADIENT Gradient for the objective function defined by
%                         the user
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        SpecifyObjectiveGradient
        
%STEPTOLERANCE Termination tolerance on the displacement in x
%
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        StepTolerance        

%TYPICALX Typical x values        
% 
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.
        TypicalX
        
%USEPARALLEL Estimate gradients in parallel
% 
%   For more information, type "doc fminimax" and see the "Options" section
%   in the FMINIMAX documentation page.        
        UseParallel
        
    end

%------------------------ Old hidden properties ---------------------------

    properties (Hidden, Dependent)

%DIAGNOSTICS Display diagnostic information 
        Diagnostics
        
%DIFFMAXCHANGE Maximum change in variables for finite-difference gradients        
        DiffMaxChange
        
%DIFFMINCHANGE Minimum change in variables for finite-difference gradients        
        DiffMinChange
        
%FINDIFFRELSTEP Scalar or vector step size factor
        FinDiffRelStep
        
%FINDIFFTYPE Finite difference type
        FinDiffType
        
%FUNVALCHECK Check whether objective function and constraints values are
%            valid
        FunValCheck
        
%GRADCONSTR Gradient for nonlinear constraint functions defined by the user
        GradConstr
        
%GRADOBJ Gradient for the objective function defined by the user
        GradObj
        
%MAXFUNEVALS Maximum number of function evaluations allowed   
        MaxFunEvals
        
%MAXITER Maximum number of iterations allowed 
        MaxIter
        
%MAXSQPITER Maximum number of SQP iterations allowed
        MaxSQPIter
        
%MERITFUNCTION Set the merit function to be used
        MeritFunction
        
%MINABSMAX Number of elements of Fi(x) to minimize the maximum absolute 
%          value of Fi
        MinAbsMax
        
%PLOTFCNS Plots various measures of progress while the algorithm executes
        PlotFcns      
        
%RELLINESRCHBND Relative bound on the line search step length
        RelLineSrchBnd
        
%RELLINESRCHBNDDURATION Number of iterations for which the bound specified 
%                       in RelLineSrchBnd should be active
        RelLineSrchBndDuration
        
%TOLCON Tolerance on the constraint violation
        TolCon
        
%TOLCONSQP Termination tolerance on inner iteration SQP constraint violation
        TolConSQP
        
%TOLFUN Termination tolerance on the function value
        TolFun

%TOLX Termination tolerance on x
        TolX
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        % This property is added in the 2nd version. We do not change the
        % Version property (base class).
        FminimaxVersion
    end
    
    % Hidden properties
    properties (Hidden, Dependent)
%NOSTOPIFFLATINFEAS If objective appears flat, only stop if feasible
%         
        NoStopIfFlatInfeas
        
%PHASEONETOTALSCALING Scale the slack variable in phase 1 of qpsub
%         
        PhaseOneTotalScaling
        
% DERIVATIVECHECK Removed from fminimax; this is to preserve backward
% compatibility
%
        DerivativeCheck        
    end
    
    properties (Hidden, Access = protected)
%OPTIONSSTORE Contains the option values and meta-data for the class
%          
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
%SOLVERNAME Name of the solver that the options are intended for
%        
        SolverName = 'fminimax';
    end
    
    properties(Hidden, Constant, GetAccess=public)
% Constant, globally visible metadata about this class.
% This data is used to spec the options in this class for internal clients
% such as: tab-complete, and the options validation
% Properties
        PropertyMetaInfo = genPropInfo();    
    end
    
    methods (Hidden)
        
        function obj = Fminimax(varargin)
%Fminimax Options for FMINIMAX
%
%   The OPTIM.OPTIONS.FMINIMAX class allows the user to create a set of
%   options for the FMINIMAX solver. For a list of options that can be set,
%   see the documentation for FMINIMAX.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX creates a set of options for FMINIMAX
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX(PARAM, VAL, ...) creates a set of options
%   for FMINIMAX with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.FMINIMAX(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.SINGLEALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS            
            
            % Call the superclass constructor
            obj = obj@optim.options.SingleAlgorithm(varargin{:});
            
            % Record the class version. Do not change Version, change
            % FminimaxVersion instead.
            obj.Version = 1;
            
            obj.FminimaxVersion = 2;
            
        end
        
        function optionFeedback = createOptionFeedback(obj)
%createOptionFeedback Create option feedback string 
%
%   optionFeedback = createOptionFeedback(obj) creates an option feedback
%   strings that are required by the extended exit messages. OPTIONFEEDBACK
%   is a structure containing strings for the options that appear in the
%   extended exit messages. These strings indicate whether the option is at
%   its 'default' value or has been 'selected'. 
            
            % It is possible for a user to pass in a vector of options to
            % the solver. Silently use the first element in this array.
            obj = obj(1);
            
            % Check if TolFunValue is the default value
            if obj.OptionsStore.SetByUser.TolFunValue
                optionFeedback.TolFunValue = 'selected';
            else
                optionFeedback.TolFunValue = 'default';
            end                
                                    
            % Check if TolFun is the default value
            if obj.OptionsStore.SetByUser.TolFun
                optionFeedback.TolFun = 'selected';
            else
                optionFeedback.TolFun = 'default';
            end
            
            % Check if TolX is the default value
            if obj.OptionsStore.SetByUser.TolX
                optionFeedback.TolX = 'selected';
            else
                optionFeedback.TolX = 'default';
            end
            
            % Check if MaxIter is the default value
            if obj.OptionsStore.SetByUser.MaxIter
                optionFeedback.MaxIter = 'selected';
            else
                optionFeedback.MaxIter = 'default';
            end
            
            % Check if MaxFunEvals is the default value
            if obj.OptionsStore.SetByUser.MaxFunEvals
                optionFeedback.MaxFunEvals = 'selected';
            else
                optionFeedback.MaxFunEvals = 'default';
            end
            
            % Check if TolCon is the default value
            if obj.OptionsStore.SetByUser.TolCon
                optionFeedback.TolCon = 'selected';
            else
                optionFeedback.TolCon = 'default';
            end
            
            % Check if MaxSQPIter is the default value
            if obj.OptionsStore.SetByUser.MaxSQPIter
                optionFeedback.MaxSQPIter = 'selected';
            else
                optionFeedback.MaxSQPIter = 'default';
            end
            
        end
        
        function obj = replaceSpecialStrings(obj)
            %replaceSpecialStrings Replace special string values 
            %
            %   obj = replaceSpecialStrings(obj) replaces special string
            %   option values with their equivalent numerical value. We
            %   currently only use this method to convert FinDiffRelStep.
            %   However, in the future we would like to move the special
            %   string replacement code from the solver files to the
            %   options classes.
           
            % Call a package function to replace string values in
            % FinDiffRelStep.
            obj.OptionsStore.Options.FinDiffRelStep = ...
                optim.options.replaceFinDiffRelStepString(obj.FinDiffRelStep);
        end
        
    end
    
    % Set/get methods
    methods

        % ---------------------- Set methods ------------------------------
        
        function obj = set.AbsoluteMaxObjectiveCount(obj, value)
            obj = setAliasProperty(obj, 'AbsoluteMaxObjectiveCount', ...
                'MinAbsMax', value);
        end

        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);        
        end         
        
        function obj = set.DerivativeCheck(obj, value)
            obj = setProperty(obj, 'DerivativeCheck', value);
        end

        function obj = set.Diagnostics(obj, value)
            obj = setProperty(obj, 'Diagnostics', value);
        end

        function obj = set.DiffMinChange(obj, value)
            obj = setProperty(obj, 'DiffMinChange', value);
        end
        
        function obj = set.DiffMaxChange(obj, value)
            obj = setProperty(obj, 'DiffMaxChange', value);
        end
        
        function obj = set.Display(obj, value)
            % Pass the possible values that the Display option can take via
            % the fourth input of setProperty.                        
            obj = setProperty(obj, 'Display', value, ...
                {'off','none','notify','notify-detailed','final', ...
                'final-detailed','iter','iter-detailed'});
        end
        
        function obj = set.FinDiffRelStep(obj, value)
            obj = setProperty(obj, 'FinDiffRelStep', value);
        end
                
        function obj = set.FiniteDifferenceStepSize(obj, value)
            obj = setAliasProperty(obj, 'FiniteDifferenceStepSize', 'FinDiffRelStep', value);
        end        

        function obj = set.FiniteDifferenceType(obj, value)
            obj = setAliasProperty(obj, 'FiniteDifferenceType', 'FinDiffType', value);     
            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            if ~obj.OptionsStore.SetByUser.FinDiffRelStep
                obj.OptionsStore.Options.FinDiffRelStep = ...
                    optim.options.getDefaultFinDiffRelStep(...
                    obj.OptionsStore.Options.FinDiffType);
            end            
        end         
        
        function obj = set.FinDiffType(obj, value)
            obj = setProperty(obj, 'FinDiffType', value);
            % If we get here, the property set has been successful and we
            % can update the OptionsStore
            if ~obj.OptionsStore.SetByUser.FinDiffRelStep
                obj.OptionsStore.Options.FinDiffRelStep = ...
                    optim.options.getDefaultFinDiffRelStep(...
                    obj.OptionsStore.Options.FinDiffType);
            end
        end

        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);        
        end         
                
        function obj = set.FunValCheck(obj, value)
            obj = setProperty(obj, 'FunValCheck', value);
        end

        function obj = set.GradConstr(obj, value)
            obj = setProperty(obj, 'GradConstr', value);
        end

        function obj = set.GradObj(obj, value)
            obj = setProperty(obj, 'GradObj', value);
        end
        
        function obj = set.MaxFunctionEvaluations(obj, value)
            obj = setAliasProperty(obj, 'MaxFunctionEvaluations', 'MaxFunEvals', value);        
        end                
        
        function obj = set.MaxFunEvals(obj, value)
            obj = setProperty(obj, 'MaxFunEvals', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);        
        end        
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end

        function obj = set.MaxSQPIter(obj, value)
            obj = setProperty(obj, 'MaxSQPIter', value);
        end
                
        function obj = set.MeritFunction(obj, value)
            obj = setProperty(obj, 'MeritFunction', value);
        end        
        
        function obj = set.MinAbsMax(obj, value)
            obj = setProperty(obj, 'MinAbsMax', value);
        end        
        
        function obj = set.NoStopIfFlatInfeas(obj, value)
            obj = setProperty(obj, 'NoStopIfFlatInfeas', value);
        end        

        function obj = set.OptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'OptimalityTolerance', 'TolFun', value);        
        end          
                
        function obj = set.OutputFcn(obj, value)
            obj = setProperty(obj, 'OutputFcn', value);
        end
        
        function obj = set.PhaseOneTotalScaling(obj, value)
            obj = setProperty(obj, 'PhaseOneTotalScaling', value);
        end

        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);        
        end          
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end
        
        function obj = set.StepTolerance(obj, value)
            obj = setAliasProperty(obj, 'StepTolerance', 'TolX', value);        
        end                
        
        function obj = set.TolX(obj, value)
            obj = setProperty(obj, 'TolX', value);
        end
        
        function obj = set.RelLineSrchBnd(obj, value)
            obj = setProperty(obj, 'RelLineSrchBnd', value);
        end
        
        function obj = set.RelLineSrchBndDuration(obj, value)
            obj = setProperty(obj, 'RelLineSrchBndDuration', value);
        end

        function obj = set.SpecifyConstraintGradient(obj, value)
            obj = setNewProperty(obj, 'SpecifyConstraintGradient', value);
        end        

        function obj = set.SpecifyObjectiveGradient(obj, value)
            obj = setNewProperty(obj, 'SpecifyObjectiveGradient', value);
        end
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end        
        
        function obj = set.TolConSQP(obj, value)
            obj = setProperty(obj, 'TolConSQP', value);
        end
        
        function obj = set.TolFun(obj, value)
            obj = setProperty(obj, 'TolFun', value);
        end
        
        function obj = set.TypicalX(obj, value)
            obj = setProperty(obj, 'TypicalX', value);
        end
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj, 'UseParallel', value);
        end

        % ---------------------- Get methods ------------------------------
        
        function value = get.AbsoluteMaxObjectiveCount(obj)
            value = obj.OptionsStore.Options.MinAbsMax;
        end

        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end                
                
        function value = get.DerivativeCheck(obj)
            value = obj.OptionsStore.Options.DerivativeCheck;
        end        
        
        function value = get.Diagnostics(obj)
            value = obj.OptionsStore.Options.Diagnostics;
        end
        
        function value = get.DiffMaxChange(obj)
            value = obj.OptionsStore.Options.DiffMaxChange;
        end
        
        function value = get.DiffMinChange(obj)
            value = obj.OptionsStore.Options.DiffMinChange;
        end
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.FiniteDifferenceStepSize(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end        
        
        function value = get.FinDiffRelStep(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end
        
        function value = get.FiniteDifferenceType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end        
        
        function value = get.FinDiffType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
                        
        function value = get.FunValCheck(obj)
            value = obj.OptionsStore.Options.FunValCheck;
        end
        
        function value = get.GradConstr(obj)
            value = obj.OptionsStore.Options.GradConstr;
        end
        
        
        function value = get.GradObj(obj)
            value = obj.OptionsStore.Options.GradObj;
        end

        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end

        function value = get.MaxFunctionEvaluations(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end        
        
        function value = get.MaxFunEvals(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end
        
        function value = get.MaxSQPIter(obj)
            value = obj.OptionsStore.Options.MaxSQPIter;
        end

        function value = get.MeritFunction(obj)
            value = obj.OptionsStore.Options.MeritFunction;
        end     
        
        function value = get.MinAbsMax(obj)
            value = obj.OptionsStore.Options.MinAbsMax;
        end     

        function value = get.NoStopIfFlatInfeas(obj)
            value = obj.OptionsStore.Options.NoStopIfFlatInfeas;
        end     

        function value = get.OptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFun;
        end        
        
        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcn;
        end
                        
        function value = get.PhaseOneTotalScaling(obj)
            value = obj.OptionsStore.Options.PhaseOneTotalScaling;
        end     
        
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end        
        
        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end
        
        function value = get.RelLineSrchBnd(obj)
            value = obj.OptionsStore.Options.RelLineSrchBnd;
        end
        
        function value = get.RelLineSrchBndDuration(obj)
            value = obj.OptionsStore.Options.RelLineSrchBndDuration;
        end

        function value = get.SpecifyConstraintGradient(obj)
            value = optim.options.OptionAliasStore.convertToLogical( ...
                obj.OptionsStore.Options.GradConstr, 'on');
        end
        
        function value = get.SpecifyObjectiveGradient(obj)
            value = optim.options.OptionAliasStore.convertToLogical( ...
                        obj.OptionsStore.Options.GradObj, 'on');
        end       

        function value = get.StepTolerance(obj)
            value = obj.OptionsStore.Options.TolX;
        end           
                
        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end
        
        function value = get.TolConSQP(obj)
            value = obj.OptionsStore.Options.TolConSQP;
        end
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFun;
        end
        
        function value = get.TolX(obj)
            value = obj.OptionsStore.Options.TolX;
        end
        
        function value = get.TypicalX(obj)
            value = obj.OptionsStore.Options.TypicalX;
        end
        
        function value = get.UseParallel(obj)
            value = obj.OptionsStore.Options.UseParallel;
        end
                
    end    
    
% Load old objects
    methods (Static = true)
        function obj = loadobj(obj)
            
            % Option DerivativeCheck was removed. It became a hidden
            % option.

            % Objects saved in R2013a will come in as structures. 
            if isstruct(obj) && obj.Version == 1
                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.Fminimax;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
   
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.FminimaxVersion = [];
                
            end
            
            % Add TolFunValue
            if isempty(obj.FminimaxVersion) || obj.FminimaxVersion < 2
                % Add TolFunValue
                obj.OptionsStore.Defaults.TolFunValue = 1e-6;
                obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                
                % Set the version number
                obj.FminimaxVersion = 2;
            end
            
            
        end
    end    
    
end

function OS = createOptionsStore
%CREATEOPTIONSSTORE Create the OptionsStore
%
%   OS = createOptionsStore creates the OptionsStore structure. This
%   structure contains the options and meta-data for option display, e.g.
%   data determining whether an option has been set by the user. This
%   function is only called when the class is first instantiated to create
%   the OptionsStore structure in its default state. Subsequent
%   instantiations of this class pick up the default OptionsStore from the
%   MCOS class definition.
%
%   Class authors must create a structure containing all the options in a
%   field of OS called Defaults. This structure must then be passed to the
%   optim.options.generateSingleAlgorithmOptionsStore function to create
%   the full OptionsStore. See below for an example for Fminimax.

% Define the option defaults for the solver
OS.Defaults.DerivativeCheck = 'off'; % hidden option
OS.Defaults.Diagnostics = 'off';
OS.Defaults.DiffMaxChange = Inf;
OS.Defaults.DiffMinChange = 0;
OS.Defaults.Display = 'final';
OS.Defaults.FinDiffRelStep = 'sqrt(eps)';
OS.Defaults.FinDiffType = 'forward';
OS.Defaults.FunValCheck = 'off';
OS.Defaults.GradConstr = 'off';
OS.Defaults.GradObj = 'off';
OS.Defaults.MaxFunEvals = '100*numberOfVariables';
OS.Defaults.MaxIter = 400;
OS.Defaults.MaxSQPIter = '10*max(numberOfVariables,numberOfInequalities+numberOfBounds)';
OS.Defaults.MeritFunction = 'multiobj';
OS.Defaults.MinAbsMax = 0;
OS.Defaults.NoStopIfFlatInfeas = 'off';
OS.Defaults.OutputFcn = [];
OS.Defaults.PhaseOneTotalScaling = 'off';
OS.Defaults.PlotFcns = [];
OS.Defaults.RelLineSrchBnd = [];
OS.Defaults.RelLineSrchBndDuration = 1;
OS.Defaults.TolCon = 1e-6;
OS.Defaults.TolConSQP = 1e-6;
OS.Defaults.TolFun = 1e-6;
OS.Defaults.TolFunValue = 1e-6;
OS.Defaults.TolX = 1e-6;
OS.Defaults.TypicalX = 'ones(numberOfVariables,1)';
OS.Defaults.UseParallel = false;

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore(OS);

end

function propInfo = genPropInfo()
% Helper function to generate constant property metadata for the Fminimax
% options class.
import optim.internal.TypeInfo;
propInfo.AbsoluteMaxObjectiveCount = TypeInfo.numericType();
propInfo.Display = TypeInfo.enumType({'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'});
propInfo.ConstraintTolerance = TypeInfo.positiveNumericType();
propInfo.FiniteDifferenceStepSize = TypeInfo.numericType();
propInfo.FiniteDifferenceType = TypeInfo.enumType({'forward','central'});
propInfo.FunctionTolerance = TypeInfo.positiveNumericType();
propInfo.MaxFunctionEvaluations = TypeInfo.integerType();
propInfo.MaxIterations = TypeInfo.integerType();
propInfo.OptimalityTolerance = TypeInfo.positiveNumericType();
propInfo.OutputFcn = TypeInfo.fcnOrEmptyType();
propInfo.PlotFcn = TypeInfo.fcnEnumType({'optimplotx', 'optimplotfunccount', 'optimplotfval','optimplotconstrviolation', 'optimplotstepsize'});
propInfo.SpecifyConstraintGradient = TypeInfo.logicalType();
propInfo.SpecifyObjectiveGradient = TypeInfo.logicalType();
propInfo.StepTolerance = TypeInfo.positiveNumericType();
propInfo.TypicalX = TypeInfo.numericType();
propInfo.UseParallel = TypeInfo.logicalType();
end
