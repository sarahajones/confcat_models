classdef (Sealed) Quadprog < optim.options.MultiAlgorithm
%

%Quadprog Options for QUADPROG
%
%   The OPTIM.OPTIONS.QUADPROG class allows the user to create a set of
%   options for the QUADPROG solver. For a list of options that can be set,
%   see the documentation for QUADPROG.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG creates a set of options for QUADPROG
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG(PARAM, VAL, ...) creates a set of options
%   for QUADPROG with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

%   Copyright 2012-2017 The MathWorks, Inc.    
    
    properties (Dependent)
%CONSTRAINTTOLERANCE Tolerance on the constraint violation
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        ConstraintTolerance    
        
%DISPLAY Level of display
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        Display

%FUNCTIONTOLERANCE Termination tolerance on the change in function value
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        FunctionTolerance        
        
%HESSIANMULTIPLYFCN Function handle for Hessian multiply function
%
%   For more information, type "doc fmincon" and see the "Options" section
%   in the QUADPROG documentation page.
        HessianMultiplyFcn        
        
%MAXITERATIONS Maximum number of iterations allowed 
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        MaxIterations

%OPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
%                    measure
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        OptimalityTolerance

%STEPTOLERANCE Termination tolerance on the change in x
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.          
        StepTolerance
        
%SUBPROBLEMALGORITHM Algorithm used to solve a subproblem
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        SubproblemAlgorithm   
                
%TYPICALX Typical x values        
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        TypicalX
        
%LINEARSOLVER Solver used by the "interior-point-convex" algorithm
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.
        LinearSolver 
    end
    
    properties (Hidden, Dependent)
        
%CONVEXITYCHECK Disable convexity checking for least squares problems
%
        ConvexCheck
        
%ENABLEPRESOLVE Enable/Disable presolve
%
        EnablePresolve
        
%DYNAMICREG Enable/Disable dynamic regularization
%
%   Applies only to dense interior-point quadprog 
        DynamicReg
        
%DIAGNOSTICS Display diagnostic information about the function
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        Diagnostics

%HESSMULT Function handle for a Hessian multiply function
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        HessMult
        
%MAXITER Maximum number of iterations allowed
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.               
        MaxIter
        
%MAXPCGITER Maximum number of PCG (preconditioned conjugate gradient)
%           iterations        
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.               
        MaxPCGIter

%PRECONDBANDWIDTH Upper bandwidth of the preconditioner for PCG    
%
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.
        PrecondBandWidth
        
%PRESOLVEOPS Operations to perform during presolve.
%
        PresolveOps
        
%TOLCON Tolerance on the constraint violation
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        TolCon
        
%TOLFUN Termination tolerance on the function value
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        TolFun
        
%TOLPCG Termination tolerance on the PCG iteration
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.        
        TolPCG
        
%TOLX Termination tolerance on x
% 
%   For more information, type "doc quadprog" and see the "Options" section
%   in the QUADPROG documentation page.          
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
        SolverName = 'quadprog';
    end
        
    properties (Hidden, SetAccess = private, GetAccess = public)
        
        % New version property added in third version
        QuadprogVersion
    end
    
    properties(Hidden, Constant, GetAccess=public)
% Constant, globally visible metadata about this class.
% This data is used to spec the options in this class for internal clients
% such as: tab-complete, and the options validation
% Properties
        PropertyMetaInfo = genPropInfo();    
    end
    
    methods (Hidden)
        
        function obj = Quadprog(varargin)
%Quadprog Options for QUADPROG
%
%   The OPTIM.OPTIONS.QUADPROG class allows the user to create a set of
%   options for the QUADPROG solver. For a list of options that can be set,
%   see the documentation for QUADPROG.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG creates a set of options for QUADPROG
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG(PARAM, VAL, ...) creates a set of options
%   for QUADPROG with the named parameters altered with the specified
%   values.
%
%   OPTS = OPTIM.OPTIONS.QUADPROG(OLDOPTS, PARAM, VAL, ...) creates a copy
%   of OLDOPTS with the named parameters altered with the specified values.
%
%   See also OPTIM.OPTIONS.MULTIALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS    

            % Call the superclass constructor
            obj = obj@optim.options.MultiAlgorithm(varargin{:});
            
            % Record the class version; Update property 'QuadprogVersion'
            % instead of superclass property 'Version'.
            obj.Version = 2;
            obj.QuadprogVersion = 7;    
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
                        
            % Check if TolCon is the default value
            if obj.OptionsStore.SetByUser.TolCon
                optionFeedback.TolCon = 'selected';
            else
                optionFeedback.TolCon = 'default';
            end
            
        end
        
    end
    
    % Set/get methods
    methods
        
        function obj = set.HessMult(obj, value)
            obj = setProperty(obj, 'HessMult', value);
        end
        
        function obj = set.HessianMultiplyFcn(obj, value)
            obj = setAliasProperty(obj, 'HessianMultiplyFcn', 'HessMult', value);
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
            obj = setNewProperty(obj, 'TolFun', value);
        end
        
        function obj = set.OptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'OptimalityTolerance', 'TolFun', value);
        end        
        
        function obj = set.FunctionTolerance(obj, value)
            obj = setAliasProperty(obj, 'FunctionTolerance', 'TolFunValue', value);
        end          
        
        function obj = set.TolCon(obj, value)
            obj = setProperty(obj, 'TolCon', value);
        end
        
        function obj = set.ConstraintTolerance(obj, value)
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);
        end        
        
        function obj = set.TypicalX(obj, value)
            obj = setProperty(obj, 'TypicalX', value);
        end
        
        function obj = set.PresolveOps(obj, value)
            obj = setProperty(obj, 'PresolveOps', value);
        end
        
        function obj = set.ConvexCheck(obj, value)
            obj = setProperty(obj, 'ConvexCheck', value);
        end
        
        function obj = set.EnablePresolve(obj, value)
            obj = setProperty(obj, 'EnablePresolve', value);
        end
        
        function obj = set.DynamicReg(obj, value)
            obj = setProperty(obj, 'DynamicReg', value);
        end
              
        function obj = set.LinearSolver(obj, value)
            obj = setNewProperty(obj, 'LinearSolver', value, obj.PropertyMetaInfo.LinearSolver.Values);
        end
        %------------------- Get functions --------------------------------
        
        function value = get.Diagnostics(obj)
            value = obj.OptionsStore.Options.Diagnostics;
        end
        
        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end
        
        function value = get.HessMult(obj)
            value = obj.OptionsStore.Options.HessMult;
        end
        
        function value = get.HessianMultiplyFcn(obj)
            value = obj.OptionsStore.Options.HessMult;
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
        
        function value = get.PrecondBandWidth(obj)
            value = obj.OptionsStore.Options.PrecondBandWidth;
        end
        
        function value = get.SubproblemAlgorithm(obj)
            value = optim.options.OptionAliasStore.mapOptionFromStore('SubproblemAlgorithm', obj.OptionsStore.Options);
        end        
        
        function obj = get.PresolveOps(obj)
            obj = obj.OptionsStore.Options.PresolveOps;
        end
        
        function obj = get.ConvexCheck(obj)
            obj = obj.OptionsStore.Options.ConvexCheck;
        end
        
        function obj = get.EnablePresolve(obj)
            obj = obj.OptionsStore.Options.EnablePresolve;
        end
        
        function obj = get.DynamicReg(obj)
            obj = obj.OptionsStore.Options.DynamicReg;
        end
        
        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end
        
        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end           
        
        function value = get.TolFun(obj)
            value = obj.OptionsStore.Options.TolFun;
        end
        
        function value = get.FunctionTolerance(obj)
            value = obj.OptionsStore.Options.TolFunValue;
        end        
        
        function value = get.OptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFun;
        end                
        
        function value = get.TolPCG(obj)
            value = obj.OptionsStore.Options.TolPCG;
        end
        
        function value = get.TolX(obj)
            value = obj.OptionsStore.Options.TolX;
        end
        
        function value = get.StepTolerance(obj)
            value = obj.OptionsStore.Options.TolX;
        end        
        
        function value = get.TypicalX(obj)
            value = obj.OptionsStore.Options.TypicalX;
        end
        
        function value = get.LinearSolver(obj)
            value = obj.OptionsStore.Options.LinearSolver;
        end
        
    end
    
    % Using mapAlgorithmName to add deprecation check
    methods (Access = protected)
        
        function name = mapAlgorithmName(obj,name)
            %MAPALGORITHMNAME Map old algorithm name to current one
            %
            %   NAME = MAPALGORITHMNAME(OBJ, OLDNAME) maps a previous name for an
            %   algorithm to its current value.
            %
            
            if strcmp(name, 'active-set')
                [linkTag,endLinkTag] = linkToAlgDefaultChangeCsh('quadprog_warn_will_error'); % links to context sensitive help
                msg = getString(message('optim:quadprog:ActSetRemoved','active-set','quadprog', ...
                    linkTag,endLinkTag,obj.OptionsStore.AlgorithmNames{:}));
                error('optim:options:QPActSetRemoved',msg);
            end
            
        end
    end    
    
    % Hidden utility methods
    methods (Hidden)
        
        function OptionsStruct = mapOptionsForSolver(~, OptionsStruct)
%mapOptionsForSolver Map structure to an optimset one
%
%   OptionsStruct = mapOptionsForSolver(obj, OptionsStruct) maps the
%   specified structure so it can be used in the solver functions and in
%   OPTIMTOOL.
%
%   This method also maps the TolFun and MaxIter options to empty if they
%   are currently marked as having a 'default dependent on problem'.

            if isfield(OptionsStruct, 'TolFun') && ...
                    strcmp(OptionsStruct.TolFun, 'default dependent on problem')
                OptionsStruct.TolFun = [];
                OptionsStruct.TolFunValue = [];
            end
            
            if isfield(OptionsStruct, 'MaxIter') && ...
                    strcmp(OptionsStruct.MaxIter, 'default dependent on problem')
                OptionsStruct.MaxIter = [];
            end    
        end

        function thisAlgorithm = createAlgorithm(obj)
%CREATEALGORITHM Create the algorithm from the options
%
%   THISALGORITHM = CREATEALGORITHM(OBJ) creates an instance of
%   optim.algorithm.QuadprogInteriorPoint from OBJ. The Options property
%   of THISALGORITHM is set to OBJ.
            switch obj.Algorithm
                case 'interior-point-convex'
                    % Create the algorithm instance and pass the options
                    % object
                    thisAlgorithm = optim.algorithm.QuadprogInteriorPoint(obj);
                otherwise
                    % For now, there aren't classes for the other quadprog
                    % algorithms
                    thisAlgorithm = [];
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
                obj = optim.options.Quadprog;
                
                % Call the superclass method to upgrade the object
                obj = upgradeFrom13a(obj, s); 
                
                % The SolverVersion property was not present in 13a. We
                % clear it here and the remainer of loadobj will set it
                % correctly.
                obj.QuadprogVersion = [];
                
            end

            % Upgrading for 13b
            if obj.Version < 2 && ...
                    ~isfield(obj.OptionsStore.AlgorithmDefaults{2}, 'PresolveOps')
                % Adding hidden property, PresolveOps
                obj.OptionsStore.AlgorithmDefaults{2}.PresolveOps = [];
                obj.OptionsStore.SetByUser.PresolveOps = false; 
                obj.OptionsStore.IsConstantDefault.PresolveOps = true; 
                obj.OptionsStore.Options.PresolveOps = [];                
            end
            
            % Upgrading for 14a
            % Change the algorithm default
            if obj.Version < 2
                % Change the default in the OptionsStore
                obj.OptionsStore.DefaultAlgorithm = 'interior-point-convex';
                
                % If the user hasn't set the Algorithm option, keep the
                % saved value of Algorithm.
                if ~obj.OptionsStore.SetByUser.Algorithm
                    obj = setPropertyNoChecks(obj, ...
                        'Algorithm', 'trust-region-reflective');
                end
                
            end
            
            % Add TolFunValue
            if isempty(obj.QuadprogVersion) || obj.QuadprogVersion < 4
                % Add TolFunValue
                obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue = 'default dependent on problem';
                obj.OptionsStore.IsConstantDefault.TolFunValue = true;
                % Set TolFunValue to whatever of TolFun was saved, but only if the selected algorithm has
                % "FunctionTolerance". Otherwise, set to its default value
                % for another algorithm
                if isfield(obj.OptionsStore.AlgorithmDefaults{obj.OptionsStore.AlgorithmIndex},'TolFunValue') && obj.OptionsStore.SetByUser.TolFun
                    obj.OptionsStore.SetByUser.TolFunValue = obj.OptionsStore.SetByUser.TolFun;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.Options.TolFun;
                else
                    obj.OptionsStore.SetByUser.TolFunValue = false;
                    obj.OptionsStore.Options.TolFunValue = obj.OptionsStore.AlgorithmDefaults{1}.TolFunValue;
                end
                
                % Objects prior to 15b are missing display-related fields
                % in OptionsStore
                obj.OptionsStore = optim.options.getDisplayOptionFieldsFor16a(...
                    obj.OptionsStore, getDefaultOptionsStore);
                
            end                       
            
            % Remove active-set
            if isempty(obj.QuadprogVersion) || obj.QuadprogVersion < 5
                % Create new options store and swap set by user values?
                if obj.OptionsStore.SetByUser.Algorithm && ...
                   strcmpi(obj.OptionsStore.Options.Algorithm,'active-set')
                    % links to context sensitive help
                    [linkTag,endLinkTag] = linkToAlgDefaultChangeCsh('quadprog_warn_will_error'); 
                    warning(message('optim:options:Quadprog:ActiveSetRemovedSwitch',...
                                    obj.OptionsStore.DefaultAlgorithm,linkTag,endLinkTag));
                    obj = obj.resetoptions('Algorithm');
                end
                newOpts = optim.options.Quadprog;
                allOptions = fieldnames(obj.OptionsStore.Options);
                for k = 1:numel(allOptions)
                    if obj.OptionsStore.SetByUser.(allOptions{k})
                        newOpts.(allOptions{k}) = obj.(allOptions{k});
                    end
                end
                obj = newOpts;
            end
            
             % Adding hidden properties, ConvexCheck, EnablePresolve and 
             % DynamicReg in 17a
            if isempty(obj.QuadprogVersion) || obj.QuadprogVersion < 6
                obj.OptionsStore.AlgorithmDefaults{2}.ConvexCheck = 'on';
                obj.OptionsStore.AlgorithmDefaults{2}.EnablePresolve = true;
                obj.OptionsStore.AlgorithmDefaults{2}.DynamicReg = 'on';

                obj.OptionsStore.SetByUser.ConvexCheck = false;
                obj.OptionsStore.SetByUser.EnablePresolve = false; 
                obj.OptionsStore.SetByUser.DynamicReg = false; 
                
                obj.OptionsStore.IsConstantDefault.ConvexCheck = true; 
                obj.OptionsStore.IsConstantDefault.EnablePresolve = true; 
                obj.OptionsStore.IsConstantDefault.DynamicReg = true; 

                obj.OptionsStore.Options.ConvexCheck = 'on';   
                obj.OptionsStore.Options.EnablePresolve = true;
                obj.OptionsStore.Options.DynamicReg = 'on';
            end
            
            % Upgrade to 18b
            if (isempty(obj.QuadprogVersion) || obj.QuadprogVersion < 7)
                % Add LinearSolver options to 'interior-point-convex' algorithm
                obj.OptionsStore.AlgorithmDefaults{2}.LinearSolver = 'auto';
            end
            
            % Set the version number
            obj.QuadprogVersion = 7;
            
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
%   the full OptionsStore. See below for an example for Quadprog.

% Define the algorithm names
OS.AlgorithmNames = {'trust-region-reflective', 'interior-point-convex'};

% Define the default algorithm
OS.DefaultAlgorithm = 'interior-point-convex';

% Define the defaults for each algorithm
% trust-region-reflective
OS.AlgorithmDefaults{1}.Diagnostics = 'off';
OS.AlgorithmDefaults{1}.Display = 'final';
OS.AlgorithmDefaults{1}.MaxIter = 'default dependent on problem';
OS.AlgorithmDefaults{1}.TolFun = 'default dependent on problem';
OS.AlgorithmDefaults{1}.TolFunValue = 'default dependent on problem';
OS.AlgorithmDefaults{1}.TolX = 100*eps;
OS.AlgorithmDefaults{1}.HessMult = [];
OS.AlgorithmDefaults{1}.MaxPCGIter = 'max(1,floor(numberOfVariables/2))';
OS.AlgorithmDefaults{1}.PrecondBandWidth = 0;
OS.AlgorithmDefaults{1}.TolPCG = 0.1;
OS.AlgorithmDefaults{1}.TypicalX = 'ones(numberOfVariables,1)';

% interior-point-convex
OS.AlgorithmDefaults{2}.Diagnostics = 'off';
OS.AlgorithmDefaults{2}.Display = 'final';
OS.AlgorithmDefaults{2}.MaxIter = 200;
OS.AlgorithmDefaults{2}.TolFun = 1e-8;
OS.AlgorithmDefaults{2}.TolX = 1e-12;
OS.AlgorithmDefaults{2}.TolCon = 1e-8;
OS.AlgorithmDefaults{2}.PresolveOps = [];
OS.AlgorithmDefaults{2}.ConvexCheck = 'on';
OS.AlgorithmDefaults{2}.EnablePresolve = true;
OS.AlgorithmDefaults{2}.DynamicReg = 'on';
OS.AlgorithmDefaults{2}.LinearSolver = 'auto';


% Call the package function to generate the OptionsStore
OS = optim.options.generateMultiAlgorithmOptionsStore(OS, 'optim.options.Quadprog');

end

function os = getDefaultOptionsStore

persistent thisos

if isempty(thisos)
    opts = optim.options.Quadprog;
    thisos = getOptionsStore(opts);
end
    
os = thisos;

end

function propInfo = genPropInfo()
% Helper function to generate constant property metadata for the Quadprog
% options class
import optim.internal.TypeInfo;
propInfo.Algorithm = TypeInfo.enumType({'interior-point-convex','trust-region-reflective'});
propInfo.ConstraintTolerance = TypeInfo.positiveNumericType();
propInfo.Display = TypeInfo.enumType({'off','final','iter','iter-detailed','final-detailed'});
propInfo.FunctionTolerance = TypeInfo.positiveNumericType();
propInfo.HessianMultiplyFcn = TypeInfo.fcnType();
propInfo.MaxIterations = TypeInfo.integerType();
propInfo.OptimalityTolerance = TypeInfo.numericType();
propInfo.StepTolerance = TypeInfo.positiveNumericType();
propInfo.SubproblemAlgorithm = TypeInfo.enumType({'cg','factorization'});
propInfo.TypicalX = TypeInfo.numericType();
propInfo.LinearSolver = TypeInfo.enumType({'auto','sparse','dense'});
end