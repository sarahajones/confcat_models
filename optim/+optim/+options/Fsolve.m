classdef (Sealed) Fsolve < optim.options.MultiAlgorithm
%

%Fsolve Options for FSOLVE
%
%   The OPTIM.OPTIONS.FSOLVE class allows the user to create a set of
%   options for the FSOLVE solver. For a list of options that can be set,
%   see the documentation for FSOLVE.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE creates a set of options for FSOLVE
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE(PARAM, VAL, ...) creates a set of options
%   for FSOLVE with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

%   Copyright 2012-2017 The MathWorks, Inc.    
    
    properties (Dependent)
        
%CHECKGRADIENTS Compare user-supplied derivatives to finite-differencing 
%               derivatives
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.       
        CheckGradients

%DISPLAY Level of display
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        Display
        
%FINITEDIFFERENCESTEPSIZE Scalar or vector step size factor
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        FiniteDifferenceStepSize
        
%FINITEDIFFERENCETYPE Finite difference type
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        FiniteDifferenceType

%FUNCTIONTOLERANCE Termination tolerance on the change in function value
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.          
        FunctionTolerance
 
%JACOBIANMULTIPLYFCN Handle to function that performs Jacobian-vector and
%                    Jacobian-matrix products
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.   
        JacobianMultiplyFcn        
        
%MAXFUNCTIONEVALUATIONS Maximum number of function evaluations allowed     
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.
        MaxFunctionEvaluations
        
%MAXITERATIONS Maximum number of iterations allowed 
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.          
        MaxIterations

%OPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
%                    measure
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.          
        OptimalityTolerance            
        
%OUTPUTFCN User-defined functions that are called at each iteration
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.         
        OutputFcn
        
%PLOTFCN Plots various measures of progress while the algorithm executes
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        PlotFcn   
        
%SPECIFYOBJECTIVEGRADIENT Specify whether a user-defined Jacobian is used
%
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.           
        SpecifyObjectiveGradient        
        
%STEPTOLERANCE Termination tolerance on x
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.          
        StepTolerance
  
%SUBPROBLEMALGORITHM Algorithm for solving a subproblem
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        SubproblemAlgorithm
        
%TYPICALX Typical x values        
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        TypicalX
        
%USEPARALLEL Option to run the finite difference calculation in parallel       
% 
%   For more information, type "doc fsolve" and see the "Options" section
%   in the FSOLVE documentation page.        
        UseParallel
    end
    
    % Hidden properties
    properties (Hidden, Dependent)
%DERIVATIVECHECK Compare user-supplied derivatives to finite-differencing 
%                derivatives
%
        DerivativeCheck
        
%DIAGNOSTICS Display diagnostic information 
%       
        Diagnostics
        
%DIFFMAXCHANGE Maximum change in variables for finite-difference gradients        
%
        DiffMaxChange
        
%DIFFMINCHANGE Minimum change in variables for finite-difference gradients        
%
        DiffMinChange
       
%FINDIFFRELSTEP Scalar or vector step size factor
%       
        FinDiffRelStep
        
%FINDIFFTYPE Finite difference type
%
        FinDiffType
        
%FUNVALCHECK Check whether objective function and constraints values are
%            valid
%
        FunValCheck

%INITDAMPING Initial Levenberg-Marquardt damping
%
        InitDamping

%JACOBIAN Specify whether a user-defined Jacobian is used
%
        Jacobian
        
%JACOBMULT Function handle for Jacobian multiply function        
%     
        JacobMult
        
%JACOBPATTERN Sparsity pattern of the Jacobian for finite differencing        
%      
        JacobPattern
               
%MAXFUNEVALS Maximum number of function evaluations allowed     
%
        MaxFunEvals
        
%MAXITER Maximum number of iterations allowed 
%
        MaxIter
        
%MAXPCGITER Maximum number of PCG (preconditioned conjugate gradient) 
%           iterations        
% 
        MaxPCGIter
              
%PLOTFCNS Plots various measures of progress while the algorithm executes
%      
        PlotFcns
        
%PRECONDBANDWIDTH Upper bandwidth of preconditioner for PCG
%   
        PrecondBandWidth
        
%SCALEPROBLEM Determine whether all constraints and the objective function
%             are normalized
% 
        ScaleProblem
        
%TOLFUN Termination tolerance on the function value
% 
        TolFun
        
%TOLPCG Termination tolerance on the PCG iteration
% 
        TolPCG
        
%TOLX Termination tolerance on x
% 
        TolX   
        
    end
    
    properties (Hidden, Access = protected)
%OPTIONSSTORE Contains the option values and meta-data for the class
%          
        OptionsStore = createOptionsStore;
    end
    
    properties (Hidden)
%SOLVERNAME Name of the solver that the options are intended for
%         
        SolverName = 'fsolve';
    end       

    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in third version
        FsolveVersion
    end
    
    properties(Hidden, Constant, GetAccess=public)
% Constant, globally visible metadata about this class.
% This data is used to spec the options in this class for internal clients
% such as: tab-complete, and the options validation
% Properties
        PropertyMetaInfo = genPropInfo();    
    end
        
    methods (Hidden)
        
        function obj = Fsolve(varargin)
%Fsolve Options for FSOLVE
%
%   The OPTIM.OPTIONS.FSOLVE class allows the user to create a set of
%   options for the FSOLVE solver. For a list of options that can be set,
%   see the documentation for FSOLVE.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE creates a set of options for FSOLVE
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE(PARAM, VAL, ...) creates a set of options
%   for FSOLVE with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.FSOLVE(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

            % Call the superclass constructor
            obj = obj@optim.options.MultiAlgorithm(varargin{:});
            
            % Record the class version; Update property 'FsolveVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.FsolveVersion = 4;    
            
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
            
            % Check if TolFun is the default value
            if obj.OptionsStore.SetByUser.TolFun
                optionFeedback.TolFun = 'selected';
            else
                optionFeedback.TolFun = 'default';
            end
            
            % Check if TolFunValue is the default value
            if obj.OptionsStore.SetByUser.TolFunValue
                optionFeedback.TolFunValue = 'selected';
            else
                optionFeedback.TolFunValue = 'default';
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

        function obj = set.ScaleProblem(obj, value)
            obj = setProperty(obj, 'ScaleProblem', value, {'none','jacobian'});
        end
        
        function obj = set.JacobPattern(obj, value)
            obj = setProperty(obj, 'JacobPattern', value);
        end
        
        function obj = set.JacobMult(obj, value)
            obj = setProperty(obj, 'JacobMult', value);
        end
        
        function obj = set.FinDiffRelStep(obj, value)
            obj = setProperty(obj, 'FinDiffRelStep', value);
        end
        
        function obj = set.FiniteDifferenceStepSize(obj, value)
            obj = setAliasProperty(obj, 'FiniteDifferenceStepSize', 'FinDiffRelStep', value);
        end        
        
        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end        
        
        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end           
        
        function obj = set.OutputFcn(obj, value)
            obj = setProperty(obj, 'OutputFcn', value);
        end        
        
        function obj = set.MaxFunEvals(obj, value)
            obj = setProperty(obj, 'MaxFunEvals', value);
        end        
        
        function obj = set.MaxFunctionEvaluations(obj, value)
            obj = setAliasProperty(obj, 'MaxFunctionEvaluations', 'MaxFunEvals', value);
        end         
        
        function obj = set.DiffMaxChange(obj, value)
            obj = setProperty(obj, 'DiffMaxChange', value);
        end        
        
        function obj = set.Jacobian(obj, value)
            obj = setProperty(obj, 'Jacobian', value);
        end        

        function obj = set.JacobianMultiplyFcn(obj, value)
            obj = setAliasProperty(obj, 'JacobianMultiplyFcn', 'JacobMult', value);
        end         
        
        function obj = set.SpecifyObjectiveGradient(obj, value)
            obj = setNewProperty(obj, 'SpecifyObjectiveGradient', value);
        end 
        
        function obj = set.DiffMinChange(obj, value)
            obj = setProperty(obj, 'DiffMinChange', value);
        end        
        
        function obj = set.DerivativeCheck(obj, value)
            obj = setProperty(obj, 'DerivativeCheck', value);
        end        
        
        function obj = set.CheckGradients(obj, value)
            obj = setNewProperty(obj, 'CheckGradients', value);
        end         
        
        function obj = set.MaxPCGIter(obj, value)
            obj = setProperty(obj, 'MaxPCGIter', value);
        end
        
        function obj = set.PrecondBandWidth(obj, value)
            obj = setProperty(obj, 'PrecondBandWidth', value);
        end
        
        function obj = set.SubproblemAlgorithm(obj, value)
            obj = setNewProperty(obj, 'SubproblemAlgorithm', value);
        end        
        
        function obj = set.TolPCG(obj, value)
            obj = setProperty(obj, 'TolPCG', value);
        end
        
        function obj = set.MaxIter(obj, value)
            obj = setProperty(obj, 'MaxIter', value);
        end
        
        function obj = set.MaxIterations(obj, value)
            obj = setAliasProperty(obj, 'MaxIterations', 'MaxIter', value);        
        end        
        
        function obj = set.TolX(obj, value)
            obj = setProperty(obj, 'TolX', value);
        end
        
        function obj = set.StepTolerance(obj, value)
            obj = setAliasProperty(obj, 'StepTolerance', 'TolX', value);        
        end                
        
        function obj = set.Display(obj, value)
            if strcmpi(value, 'testing')
                % Set Display to the undocumented value, 'testing'.
                obj = setPropertyNoChecks(obj, 'Display', 'testing');                
            else
                % Pass the possible values that the Display option can take
                % via the fourth input of setProperty.
                obj = setProperty(obj, 'Display', value, ...
                    {'off','none','final', ...
                    'final-detailed','iter','iter-detailed'});
            end
        end
        
        function obj = set.Diagnostics(obj, value)
            obj = setProperty(obj, 'Diagnostics', value);
        end
        
        function obj = set.TolFun(obj, value)
            % TolFun maps to FunctionTolerance for levenberg-marquardt and
            % to OptimalityTolerance for every other algorithm
            if strcmp(obj.Algorithm, 'levenberg-marquardt')
                obj = setNewProperty(obj, 'TolFunValue', value);
            else
                obj = setNewProperty(obj, 'TolFun', value);
            end
        end
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);        
        end  
                
        function obj = set.OptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'OptimalityTolerance', 'TolFun', value);        
        end         
        
        function obj = set.TypicalX(obj, value)
            obj = setProperty(obj, 'TypicalX', value);
        end
        
        function obj = set.UseParallel(obj, value)
            obj = setProperty(obj,'UseParallel', value);
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
        
        function obj = set.FunValCheck(obj, value)
            obj = setProperty(obj, 'FunValCheck', value);
        end
        
        function obj = set.InitDamping(obj, value)
            obj = setProperty(obj, 'InitDamping', value);
        end
             
        %---------------- Get functions -----------------------------------
        
        function value = get.DerivativeCheck(obj)
            value = obj.OptionsStore.Options.DerivativeCheck;
        end
        
        function value = get.CheckGradients(obj)
            value = optim.options.OptionAliasStore.convertToLogical( ...
                        obj.OptionsStore.Options.DerivativeCheck, 'on');
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

        function value = get.FinDiffType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end
        
        function value = get.FiniteDifferenceType(obj)
            value = obj.OptionsStore.Options.FinDiffType;
        end        
        
        function value = get.FinDiffRelStep(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end
        
        function value = get.FiniteDifferenceStepSize(obj)
            value = obj.OptionsStore.Options.FinDiffRelStep;
        end        
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end         
        
        function value = get.FunValCheck(obj)
            value = obj.OptionsStore.Options.FunValCheck;
        end
        
        function value = get.InitDamping(obj)
            value = obj.OptionsStore.Options.InitDamping;
        end  
        
        function value = get.Jacobian(obj)
            value = obj.OptionsStore.Options.Jacobian;
        end

        function value = get.JacobianMultiplyFcn(obj)
            value = obj.OptionsStore.Options.JacobMult;
        end         
        
        function value = get.JacobMult(obj)
            value = obj.OptionsStore.Options.JacobMult;
        end        

        function value = get.JacobPattern(obj)
            value = obj.OptionsStore.Options.JacobPattern;
        end        
         
        function value = get.MaxFunEvals(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end
        
        function value = get.MaxFunctionEvaluations(obj)
            value = obj.OptionsStore.Options.MaxFunEvals;
        end                 
        
        function value = get.MaxIter(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end
        
        function value = get.MaxIterations(obj)
            value = obj.OptionsStore.Options.MaxIter;
        end        
        
        function value = get.MaxPCGIter(obj)
            value = obj.OptionsStore.Options.MaxPCGIter;
        end
        
        function value = get.OptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFun;
        end            
        
        function value = get.OutputFcn(obj)
            value = obj.OptionsStore.Options.OutputFcn;
        end
        
        function value = get.PlotFcn(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end          
        
        function value = get.PlotFcns(obj)
            value = obj.OptionsStore.Options.PlotFcns;
        end         
        
        function value = get.PrecondBandWidth(obj)
            value = obj.OptionsStore.Options.PrecondBandWidth;
        end       

        function value = get.ScaleProblem(obj)
            value = obj.OptionsStore.Options.ScaleProblem;
        end        
        
        function value = get.SpecifyObjectiveGradient(obj)
            value = optim.options.OptionAliasStore.convertToLogical( ...
                        obj.OptionsStore.Options.Jacobian, 'on');
        end           
        
        function value = get.StepTolerance(obj)
            value = obj.OptionsStore.Options.TolX;
        end           
        
        function value = get.SubproblemAlgorithm(obj)
            value = optim.options.OptionAliasStore.mapOptionFromStore('SubproblemAlgorithm', obj.OptionsStore.Options);
        end        
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFun;
        end
        
        function value = get.TolPCG(obj)
            value = obj.OptionsStore.Options.TolPCG;
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

    % Map 'trust-region-reflective' to 'trust-region'
    methods (Access = protected)
        
        function name = mapAlgorithmName(~, name)
            %MAPALGORITHMNAME Map old algorithm name to current one
            %
            %   NAME = MAPALGORITHMNAME(OBJ, OLDNAME) maps a previous name for an
            %   algorithm to its current value.
            %
            
            if strcmp(name, 'trust-region-reflective')
                name = 'trust-region';
            end
            
        end
    end
    
    
    % Load old objects
    methods (Static = true)
        function obj = loadobj(obj)
            
            % Objects saved in R2013a will come in as structures. 
            if isstruct(obj) && obj.Version == 1

                % Save the existing structure
                s = obj;
                
                % Create a new object
                obj = optim.options.Fsolve;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.FsolveVersion = [];
                
            end        
            
            % Add the Levenberg-Marquardt parameter. Also, update
            % FinDiffRelStep default that was upgraded in R2013b.            
            if obj.Version < 2
                
                % Add the InitDamping option to the OptionsStore
                obj.OptionsStore.AlgorithmDefaults{3}.InitDamping = 0.01; 
                obj.OptionsStore.SetByUser.InitDamping = false; 
                obj.OptionsStore.IsConstantDefault.InitDamping = true; 
                obj.OptionsStore.Options.InitDamping = 0.01;
                
                % Upgrade FinDiffRelStep default
                obj.OptionsStore.AlgorithmDefaults{1}.FinDiffRelStep = 'sqrt(eps)';
                obj.OptionsStore.AlgorithmDefaults{2}.FinDiffRelStep = 'sqrt(eps)';
                obj.OptionsStore.AlgorithmDefaults{3}.FinDiffRelStep = 'sqrt(eps)';
                if ~isSetByUser(obj, 'FinDiffRelStep')
                    obj.OptionsStore.Options.FinDiffRelStep = 'sqrt(eps)';
                end     
                
            end

            % Set the version number
            obj.Version = 2;
            
            % Upgrading to 15b
            % Introduce FsolveVersion field
            if isempty(obj.FsolveVersion)
                % Use 'FsolveVersion' property instead of 'Version' property
                % because 'Version' is for the superclass and 'FsolveVersion' is
                % for this (derived) class. However, 'FsolveVersion' was added in
                % the third version, we check only for the second version and
                % add this property. For all other version, check only the
                % 'FsolveVersion' property.
                obj.FsolveVersion = 2; % update object
            end
            
            if obj.FsolveVersion < 3
                % Set the version number
                obj.FsolveVersion = 3;
            end
            
            % Upgrading to 16a
            if obj.FsolveVersion < 4
                
                % Add the UseParallel option to the OptionsStore
                obj.OptionsStore.AlgorithmDefaults{1}.UseParallel = false;
                obj.OptionsStore.AlgorithmDefaults{2}.UseParallel = false;
                obj.OptionsStore.AlgorithmDefaults{3}.UseParallel = false;
                obj.OptionsStore.SetByUser.UseParallel = false;
                obj.OptionsStore.IsConstantDefault.UseParallel = true;
                obj.OptionsStore.Options.UseParallel = false;
                
                % Add TolFunValue
                obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue = 1e-6;
                obj.OptionsStore.AlgorithmDefaults{2}.TolFunValue = 1e-6;
                obj.OptionsStore.AlgorithmDefaults{3}.TolFunValue = 1e-6;
                obj.OptionsStore.IsConstantDefault.TolFunValue = true;
                
                % Set TolFunValue to whatever of TolFun was saved, but only
                % if the selected algorithm has "FunctionTolerance".
                % Otherwise, set to its default value for another algorithm
                if isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex},'TolFunValue') && obj.OptionsStore.SetByUser.TolFun
                    obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                else
                    obj.OptionsStore.SetByUser.TolFunValue = false;
                    obj.OptionsStore.Options.TolFunValue = 1e-6;
                end
                
                % Remove TolFun field for Levenberg-Marquardt
                obj.OptionsStore.AlgorithmDefaults{3} = ...
                    rmfield(obj.OptionsStore.AlgorithmDefaults{3}, 'TolFun');
                
                % Objects prior to 15b are missing display-related fields
                % in OptionsStore
                obj.OptionsStore = optim.options.getDisplayOptionFieldsFor16a(...
                    obj.OptionsStore, getDefaultOptionsStore);
                
                % Renaming 'trust-region-reflective' to 'trust-region'
                obj.OptionsStore.AlgorithmNames{2} = 'trust-region';
                if strcmp(obj.OptionsStore.Options.Algorithm, 'trust-region-reflective')
                    obj.OptionsStore.Options.Algorithm = 'trust-region';                    
                end
                
                % Set the version number
                obj.FsolveVersion = 4;
                
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
%   Class authors must create a structure containing the following fields:-
%
%   AlgorithmNames   : Cell array of algorithm names for the solver
%   DefaultAlgorithm : String containing the name of the default algorithm
%   AlgorithmDefaults: Cell array of structures. AlgorithmDefaults{i}
%                      holds a structure containing the defaults for 
%                      AlgorithmNames{i}.
%
%   This structure must then be passed to the
%   optim.options.generateMultiAlgorithmOptionsStore function to create
%   the full OptionsStore. See below for an example for Fsolve.

% Define the algorithm names
OS.AlgorithmNames = {'trust-region-dogleg', 'trust-region', 'levenberg-marquardt'};

% Define the default algorithm
OS.DefaultAlgorithm = 'trust-region-dogleg';

% Define the defaults for each algorithm
% trust-region-dogleg
OS.AlgorithmDefaults{1}.DerivativeCheck = 'off';
OS.AlgorithmDefaults{1}.Diagnostics = 'off';
OS.AlgorithmDefaults{1}.DiffMaxChange = Inf;
OS.AlgorithmDefaults{1}.DiffMinChange = 0;
OS.AlgorithmDefaults{1}.Display = 'final';
OS.AlgorithmDefaults{1}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{1}.FinDiffType = 'forward';
OS.AlgorithmDefaults{1}.FunValCheck = 'off';
OS.AlgorithmDefaults{1}.Jacobian = 'off';
OS.AlgorithmDefaults{1}.MaxFunEvals = '100*numberOfVariables';
OS.AlgorithmDefaults{1}.MaxIter = 400;
OS.AlgorithmDefaults{1}.OutputFcn = [];
OS.AlgorithmDefaults{1}.PlotFcns = [];
OS.AlgorithmDefaults{1}.TolFun = 1e-6;
OS.AlgorithmDefaults{1}.TolFunValue = 1e-6;
OS.AlgorithmDefaults{1}.TolX = 1e-6;
OS.AlgorithmDefaults{1}.TypicalX = 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{1}.UseParallel = false;

% trust-region
OS.AlgorithmDefaults{2}.DerivativeCheck = 'off';
OS.AlgorithmDefaults{2}.Diagnostics = 'off';
OS.AlgorithmDefaults{2}.DiffMaxChange = Inf;
OS.AlgorithmDefaults{2}.DiffMinChange = 0;
OS.AlgorithmDefaults{2}.Display = 'final';
OS.AlgorithmDefaults{2}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{2}.FinDiffType = 'forward';
OS.AlgorithmDefaults{2}.FunValCheck = 'off';
OS.AlgorithmDefaults{2}.Jacobian = 'off';
OS.AlgorithmDefaults{2}.MaxFunEvals = '100*numberOfVariables';
OS.AlgorithmDefaults{2}.MaxIter = 400;
OS.AlgorithmDefaults{2}.OutputFcn = [];
OS.AlgorithmDefaults{2}.PlotFcns = [];
OS.AlgorithmDefaults{2}.TolFun = 1e-6;
OS.AlgorithmDefaults{2}.TolFunValue = 1e-6;
OS.AlgorithmDefaults{2}.TolX = 1e-6;
OS.AlgorithmDefaults{2}.TypicalX = 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{2}.UseParallel = false;
OS.AlgorithmDefaults{2}.JacobMult = [];
OS.AlgorithmDefaults{2}.JacobPattern = 'sparse(ones(Jrows,Jcols))';
OS.AlgorithmDefaults{2}.MaxPCGIter = 'max(1,floor(numberOfVariables/2))';
OS.AlgorithmDefaults{2}.PrecondBandWidth = Inf;
OS.AlgorithmDefaults{2}.TolPCG = 0.1;

% levenberg-marquardt
OS.AlgorithmDefaults{3}.DerivativeCheck = 'off';
OS.AlgorithmDefaults{3}.Diagnostics = 'off';
OS.AlgorithmDefaults{3}.DiffMaxChange = Inf;
OS.AlgorithmDefaults{3}.DiffMinChange = 0;
OS.AlgorithmDefaults{3}.Display = 'final';
OS.AlgorithmDefaults{3}.FinDiffRelStep = 'sqrt(eps)';
OS.AlgorithmDefaults{3}.FinDiffType = 'forward';
OS.AlgorithmDefaults{3}.FunValCheck = 'off';
OS.AlgorithmDefaults{3}.Jacobian = 'off';
OS.AlgorithmDefaults{3}.InitDamping = 0.01;
OS.AlgorithmDefaults{3}.MaxFunEvals = '200*numberOfVariables';
OS.AlgorithmDefaults{3}.MaxIter = 400;
OS.AlgorithmDefaults{3}.OutputFcn = [];
OS.AlgorithmDefaults{3}.PlotFcns = [];
OS.AlgorithmDefaults{3}.TolFunValue = 1e-6;
OS.AlgorithmDefaults{3}.TolX = 1e-6;
OS.AlgorithmDefaults{3}.TypicalX = 'ones(numberOfVariables,1)';
OS.AlgorithmDefaults{3}.UseParallel = false;
OS.AlgorithmDefaults{3}.ScaleProblem = 'none';

% Call the package function to generate the OptionsStore
OS = optim.options.generateMultiAlgorithmOptionsStore(OS, 'optim.options.Fsolve');

end

function os = getDefaultOptionsStore

persistent thisos

if isempty(thisos)
    opts = optim.options.Fsolve;
    thisos = getOptionsStore(opts);
end
    
os = thisos;

end

function propInfo = genPropInfo()
% Helper function to generate constant property metadata for the Fsolve
% options class.
import optim.internal.TypeInfo;

propInfo.Algorithm = TypeInfo.enumType({'trust-region-dogleg','trust-region','levenberg-marquardt'});
propInfo.CheckGradients = TypeInfo.logicalType();
propInfo.Display = TypeInfo.enumType({'off','iter','iter-detailed','final','final-detailed'});
propInfo.FiniteDifferenceStepSize = TypeInfo.numericType();
propInfo.FiniteDifferenceType = TypeInfo.enumType({'forward','central'});
propInfo.FunctionTolerance = TypeInfo.positiveNumericType();
propInfo.JacobianMultiplyFcn = TypeInfo.fcnOrEmptyType();
propInfo.MaxFunctionEvaluations = TypeInfo.integerType();
propInfo.MaxIterations = TypeInfo.integerType();
propInfo.OptimalityTolerance = TypeInfo.positiveNumericType();
propInfo.OutputFcn = TypeInfo.fcnType();
propInfo.PlotFcn = TypeInfo.fcnEnumType({'optimplotx', 'optimplotfunccount', 'optimplotfval', 'optimplotstepsize', 'optimplotfirstorderopt'});
propInfo.SpecifyObjectiveGradient = TypeInfo.logicalType();
propInfo.StepTolerance = TypeInfo.positiveNumericType();
propInfo.SubproblemAlgorithm = TypeInfo.enumType({'cg','factorization'});
propInfo.TypicalX = TypeInfo.numericType();
propInfo.UseParallel = TypeInfo.logicalType();
end