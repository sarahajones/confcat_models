classdef (Sealed) Intlinprog < optim.options.SingleAlgorithm
    
%Intlinprog Options for INTLINPROG
%
%   The OPTIM.OPTIONS.INTLINPROG class allows the user to create a set of
%   options for the INTLINPROG solver. For a list of options that can be
%   set, see the documentation for INTLINPROG.
%
%   OPTS = OPTIM.OPTIONS.INTLINPROG creates a set of options for INTLINPROG
%   with the options set to their default values.
%
%   OPTS = OPTIM.OPTIONS.INTLINPROG(PARAM, VAL, ...) creates a set of
%   options for INTLINPROG with the named parameters altered with the
%   specified values.
%
%   OPTS = OPTIM.OPTIONS.INTLINPROG(OLDOPTS, PARAM, VAL, ...) creates a
%   copy of OLDOPTS with the named parameters altered with the specified
%   values.
%
%   See also OPTIM.OPTIONS.SINGLEALGORITHM, OPTIM.OPTIONS.SOLVEROPTIONS

    %   Copyright 2013-2017 The MathWorks, Inc.
    
    properties (Dependent)
                
%ABSOLUTEGAPTOLERANCE Maximum absolute gap between the lower and upper
%bound estimates of the solution (a non-negative scalar). The default is 0.
        AbsoluteGapTolerance
        
%BRANCRULE Rule the algorithm uses to select the branch variable in the
%search tree. The choices are 'maxfun', 'mostfractional', 'maxpscost' (default)
        BranchRule
        
%CONSTRAINTTOLERANCE Tolerance on the constraint violation (a positive
%scalar). The default is 1e-4.
        ConstraintTolerance
        
%CUTGENERATION Method the algorithm uses for generating extra linear
%constraints in the problem. The choices are 'none', 'basic',
%'intermediate', or 'advanced'. A more advanced choice of CutGeneration
%indicates that more constraints will be added to the problem and/or more
%elaborate cut generation algorithms will be used.
        CutGeneration
        
%CUTMAXITERATIONS Maximum number of passes of cut generation performed. The
%default is 10.
        CutMaxIterations        

%DISPLAY Level of display. The choices are 'off', 'iter and 'final'
%(default).
        Display

%HEURISTICS Specify how the algorithm generates upper bound estimates of
%the objective. The choices are 'none', 'basic', 'intermediate', 'advanced'.
        Heuristics

%HEURISTICSMAXNODES Maximum number of nodes that the heuristics will be
%applied for (a positive integer). The default is 50.
        HeuristicsMaxNodes

%INTEGERPREPROCESS Integer preprocessing. The choices are 'none', 'basic'
%(default) or 'advanced'.
        IntegerPreprocess

%INTEGERTOLERANCE Tolerance within which the value of a variable is considered
%to be integral (a positive scalar). The default is 1e-5.
        IntegerTolerance

%LPMAXITERATIONS Maximum number of iterations the LP-solver performs to
%solve the LP-relaxation problem at each node (a positive integer). The
%default is 'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))'.
        LPMaxIterations

%LPOPTIMALITYTOLERANCE Termination tolerance on the first-order optimality
%measure of a linear programming relaxation problem (a positive scalar).
%The default is 1e-7.
        LPOptimalityTolerance

%MAXNODES Maximum number of nodes, the function searches (a positive
%integer). The default is 1e7.
        MaxNodes

%MAXFEASIBLEPOINTS Maximum number of integer solutions that intlinprog will
%find (a positive integer). The default is Inf.
        MaxFeasiblePoints

%MAXTIME Maximum amount of CPU time in seconds the function runs (a
%positive scalar). The default value is 7200.
        MaxTime

%NODESELECTION Strategy the algorithm uses to select the next node to
%search in the search tree. The choices are 'minobj', 'mininfeas',
%'simplebestproj' (default).
        NodeSelection

%OBJECTIVECUTOFF A node is not analyzed if the upper bound for the node is
%above this value (scalar). The default is Inf.
        ObjectiveCutOff

%OBJECTIVEIMPROVEMENTTHRESHOLD Relative improvement required to accept a
%new integer feasible solution (a positive scalar). The default is 0.0.
        ObjectiveImprovementThreshold

%OUTPUTFCN A function handle or a cell array of function handles. The
%default is [].
        OutputFcn

%PLOTFCN A function handle or a cell array of function handles. The
%default is [].
        PlotFcn

%RELATIVEGAPTOLERANCE Maximum relative gap between the lower and upper bound estimates
%of the solution (a non-negative scalar). The default is 1e-4
        RelativeGapTolerance

%ROOTLPALGORITHM Algorithm used to solve the problem at the root node.
%The choices are 'primal-simplex' or 'dual-simplex' (default).
        RootLPAlgorithm

%ROOTLPMAXITERATIONS Maximum number of iterations that the root node
%algorithm can take to solve the root node problem (a positive integer).
%The default is 'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))'.
        RootLPMaxIterations
    end

    % Hidden properties
    properties(Hidden, Dependent)

%BRANCHINGRULE Rule the algorithm uses to select the branch variable in the
%search tree. The choices are  'maxfun', 'mostfractional', 'maxpscost' (default)
%'strongpscost' and 'reliability'.
        BranchingRule

%CUTGENMAXITER Maximum number of passes of cut generation performed.
%The default is 10.
        CutGenMaxIter

%IPPREPROCESS Integer preprocessing. The choices are 'none', 'normal'
%(default) or 'advanced'.
        IPPreprocess

%LPMAXITER Maximum number of iterations the LP-solver performs to solve the
%LP-relaxation problem at each node (a positive integer). The default is
%'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))'.
        LPMaxIter

%LPPREPROCESS Preprocessing for the linear problems. The choices are
%'none', or 'basic' (default).
        LPPreprocess

%MAXNUMFEASPOINTS Maximum number of integer solutions that intlinprog will
%find (a positive integer). The default is Inf.
        MaxNumFeasPoints

%PLOTFCNS A function handle or a cell array of function handles. The
%default is [].
        PlotFcns

%RELOBJTHRESHOLD Relative improvement required to accept a new integer
%feasible solution (a positive scalar). The default is 0.0.
        RelObjThreshold

%ROOTLPMAXITER Maximum number of iterations that the root node algorithm
%can take to solve the root node problem (a positive integer). The default
%is 'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))'.
        RootLPMaxIter

%TOLCON Tolerance on the constraint violation (a positive scalar). The
%default is 1e-4.
        TolCon

%TOLFUNLP Termination tolerance on the function value of a linear
%programming relaxation problem (a positive scalar). The default is 1e-4.
        TolFunLP

%TOLGAPABS Maximum absolute gap between the lower and upper bound estimates
%of the solution (a non-negative scalar). The default is 0.
        TolGapAbs

%TOLGAPREL Maximum relative gap between the lower and upper bound estimates
%of the solution (a non-negative scalar). The default is 1e-4
        TolGapRel

%TOLINTEGER Tolerance within which the value of a variable is considered
%to be integral (a positive scalar). The default is 1e-5.
        TolInteger
    end

    properties (SetAccess = private, GetAccess = private)

        %INTERNALOPTIONS Internal options
        InternalOptions = struct();
        % This property is added in 2nd version. We do not change the
        % Version property (base class).
        IntlinprogVersion
    end


    properties (Hidden, Access = protected)

        %OPTIONSSTORE Contains the option values and meta-data for the class
        %
        OptionsStore = createOptionsStore;

    end

    properties (Hidden)

        %SOLVERNAME Name of the solver that the options are intended for
        %
        SolverName = 'intlinprog';

    end
    
    properties(Hidden, Constant, GetAccess=public)
% Constant, globally visible metadata about this class.
% This data is used to spec the options in this class for internal clients
% such as: tab-complete, and the options validation
% Properties
        PropertyMetaInfo = genPropInfo();    
    end

    methods (Hidden)

        function obj = Intlinprog(varargin)

            % Call the superclass constructor
            obj = obj@optim.options.SingleAlgorithm(varargin{:});

            % Record the class version. Do not change Version, change
            % IntlinprogVersion instead.
            obj.Version = 1;
            obj.IntlinprogVersion = 5;
        end

    end

    methods (Static = true)

      function obj = loadobj(obj)
        if isempty(obj.IntlinprogVersion)
          % Use IntlinprogVersion property instead of Version property
          % because Version is for the super class only. However,
          % IntlinprogVersion was added in 2nd version of this class. We
          % check only for the 1st version and add this property. For all
          % other version, check only IntlinprogVersion property.
          obj.IntlinprogVersion = 1; % update object
        end

        if obj.IntlinprogVersion < 2
          % OutputFcn was added for IntlinprogVersion == 2
           obj.OptionsStore.Defaults.OutputFcn = [];
           obj.OptionsStore.SetByUser.OutputFcn = false;
           obj.OptionsStore.Options.OutputFcn = [];

           % PlotFcns was added for IntlinprogVersion == 2
           obj.OptionsStore.Defaults.PlotFcns = [];
           obj.OptionsStore.SetByUser.PlotFcns = false;
           obj.OptionsStore.Options.PlotFcns = [];
        end

        % Upgrading to 17a
        if obj.IntlinprogVersion < 4
                % Change the default in the OptionsStore
                obj.OptionsStore.Defaults.RootLPMaxIter = ...
                    'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))';

                % If the user hasn't set the RootLPMaxIterations option, keep the
                % saved value of RootLPMaxIterations.
                if ~obj.OptionsStore.SetByUser.RootLPMaxIter
                    obj = setPropertyNoChecks(obj, ...
                        'RootLPMaxIterations', 30000);
                end
        end

        % Upgrading to 18a
        if obj.IntlinprogVersion < 5
                % Change the default in the OptionsStore
                obj.OptionsStore.Defaults.RelObjThreshold = 0;

                % If the user hasn't set the ObjectiveImprovementThreshold
                % option, keep the saved value of
                % ObjectiveImprovementThreshold.
                if ~obj.OptionsStore.SetByUser.RelObjThreshold
                    obj = setPropertyNoChecks(obj, ...
                        'ObjectiveImprovementThreshold', 1e-4);
                end
                
                % Change the default in the OptionsStore
                obj.OptionsStore.Defaults.LPMaxIter = ...
                    'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))';

                % If the user hasn't set the LPMaxIterations option, keep the
                % saved value of LPMaxIterations.
                if ~obj.OptionsStore.SetByUser.LPMaxIter
                    obj = setPropertyNoChecks(obj, ...
                        'LPMaxIterations', 30000);
                end                
        end

        % Update the version info to the current number
        obj.IntlinprogVersion = 5;
      end

    end


    % Set methods
    methods

        function obj = set.BranchingRule(obj, value)
            obj = setProperty(obj, 'BranchingRule', value, ...
                {'mostfractional', 'maxpscost', 'maxfun', 'strongpscost', 'reliability'});
        end

        function obj = set.BranchRule(obj, value)
            obj = setAliasProperty(obj, 'BranchRule', 'BranchingRule', value, ...
                {'mostfractional', 'maxpscost', 'maxfun', 'strongpscost', 'reliability'});
        end

        function obj = set.CutGeneration(obj, value)
            obj = setProperty(obj, 'CutGeneration', value, ...
                {'none', 'basic', 'intermediate', 'advanced'});
        end

        function obj = set.CutGenMaxIter(obj, value)
            obj = setProperty(obj, 'CutGenMaxIter', value);
        end

        function obj = set.CutMaxIterations(obj, value)
            obj = setAliasProperty(obj, 'CutMaxIterations', ...
                                   'CutGenMaxIter', value);
        end

        function obj = set.Display(obj, value)
            % Pass the possible values that the Display option can take via
            % the fourth input of setProperty.
            obj = setProperty(obj, 'Display', value, ...
                {'off', 'none', 'final', 'iter'});
        end

        function obj = set.HeuristicsMaxNodes(obj, value)
            obj = setProperty(obj, 'HeuristicsMaxNodes', value);
        end

        function obj = set.Heuristics(obj, value)
            obj = setProperty(obj, 'Heuristics', value, ...
                {'none', 'basic', 'intermediate', 'advanced','rins','rss', ...
                'round','diving','rss-diving','rins-diving','round-diving'});
        end

        function obj = set.IPPreprocess(obj, value)
            obj = setProperty(obj, 'IPPreprocess', value, ...
                {'none', 'basic', 'advanced'});
        end

        function obj = set.IntegerPreprocess(obj, value)
            obj = setAliasProperty(obj, 'IntegerPreprocess', 'IPPreprocess', value, ...
                {'none', 'basic', 'advanced'});
        end

        function obj = set.LPMaxIter(obj, value)
            obj = setProperty(obj, 'LPMaxIter', value);
        end

        function obj = set.LPMaxIterations(obj, value)
            obj = setAliasProperty(obj, 'LPMaxIterations', 'LPMaxIter', value);
        end

        function obj = set.LPPreprocess(obj, value)
            obj = setProperty(obj, 'LPPreprocess', value, ...
                {'none', 'basic'});
        end

        function obj = set.MaxNodes(obj, value)
            obj = setProperty(obj, 'MaxNodes', value);
        end

        function obj = set.MaxNumFeasPoints(obj, value)
            obj = setProperty(obj, 'MaxNumFeasPoints', value);
        end

        function obj = set.MaxFeasiblePoints(obj, value)
            obj = setAliasProperty(obj, 'MaxFeasiblePoints', 'MaxNumFeasPoints', value);
        end

        function obj = set.MaxTime(obj, value)
            obj = setProperty(obj, 'MaxTime', value);
        end

        function obj = set.NodeSelection(obj, value)
            obj = setProperty(obj, 'NodeSelection', value, ...
                {'minobj', 'mininfeas', 'simplebestproj'});
        end

        function obj = set.ObjectiveCutOff(obj, value)
            obj = setProperty(obj, 'ObjectiveCutOff', value);
        end

        function obj = set.OutputFcn(obj, value)
            obj = setProperty(obj, 'OutputFcn', value);
        end

        function obj = set.PlotFcns(obj, value)
            obj = setProperty(obj, 'PlotFcns', value);
        end

        function obj = set.PlotFcn(obj, value)
            obj = setAliasProperty(obj, 'PlotFcn', 'PlotFcns', value);
        end

        function obj = set.RootLPAlgorithm(obj, value)
            obj = setProperty(obj, 'RootLPAlgorithm', value, ...
                 {'primal-simplex', 'dual-simplex'});
        end

        function obj = set.RootLPMaxIter(obj, value)
            obj = setProperty(obj, 'RootLPMaxIter', value);
        end

        function obj = set.RootLPMaxIterations(obj, value)
            obj = setAliasProperty(obj, 'RootLPMaxIterations', 'RootLPMaxIter', value);
        end

        function obj = set.TolCon(obj, value)
            % We have to perform this check here, because the TolCon check
            % in optim.options.checkfield only ensures that TolCon is a
            % non-negative real.
            if value < 1e-9 || value > 1e-3
                error(message('optim:options:Intlinprog:TolConOutOfRange','TolCon'));
            end
            obj = setProperty(obj, 'TolCon', value);
        end

        function obj = set.ConstraintTolerance(obj, value)
            % We have to perform this check here, because the TolCon check
            % in optim.options.checkfield only ensures that TolCon is a
            % non-negative real.
            if value < 1e-9 || value > 1e-3
                error(message('optim:options:Intlinprog:TolConOutOfRange','ConstraintTolerance'));
            end
            obj = setAliasProperty(obj, 'ConstraintTolerance', 'TolCon', value);
        end

        function obj = set.RelObjThreshold(obj, value)
            obj = setProperty(obj, 'RelObjThreshold', value);
        end

        function obj = set.ObjectiveImprovementThreshold(obj, value)
            obj = setAliasProperty(obj, 'ObjectiveImprovementThreshold', 'RelObjThreshold', value);
        end

        function obj = set.TolFunLP(obj, value)
            obj = setProperty(obj, 'TolFunLP', value);
        end

        function obj = set.LPOptimalityTolerance(obj, value)
            obj = setAliasProperty(obj, 'LPOptimalityTolerance', 'TolFunLP', value);
        end

        function obj = set.TolGapAbs(obj, value)
            obj = setProperty(obj, 'TolGapAbs', value);
        end

        function obj = set.AbsoluteGapTolerance(obj, value)
            obj = setAliasProperty(obj, 'AbsoluteGapTolerance', 'TolGapAbs', value);
        end

        function obj = set.TolGapRel(obj, value)
            obj = setProperty(obj, 'TolGapRel', value);
        end

        function obj = set.RelativeGapTolerance(obj, value)
            obj = setAliasProperty(obj, 'RelativeGapTolerance', 'TolGapRel', value);
        end

        function obj = set.TolInteger(obj, value)
            obj = setProperty(obj, 'TolInteger', value);
        end

        function obj = set.IntegerTolerance(obj, value)
            obj = setAliasProperty(obj, 'IntegerTolerance', 'TolInteger', value);
        end
    end
    % Get methods
    methods

        function value = get.BranchRule(obj)
            value = obj.OptionsStore.Options.BranchingRule;
        end

        function value = get.BranchingRule(obj)
            value = obj.OptionsStore.Options.BranchingRule;
        end

        function value = get.CutGeneration(obj)
            value = obj.OptionsStore.Options.CutGeneration;
        end

        function value = get.CutGenMaxIter(obj)
            value = obj.OptionsStore.Options.CutGenMaxIter;
        end

        function value = get.CutMaxIterations(obj)
            value = obj.OptionsStore.Options.CutGenMaxIter;
        end

        function value = get.Display(obj)
            value = obj.OptionsStore.Options.Display;
        end

        function value = get.Heuristics(obj)
            value =  obj.OptionsStore.Options.Heuristics;
        end

        function value = get.HeuristicsMaxNodes(obj)
            value = obj.OptionsStore.Options.HeuristicsMaxNodes;
        end

        function value = get.IPPreprocess(obj)
            value = obj.OptionsStore.Options.IPPreprocess;
        end

        function value = get.IntegerPreprocess(obj)
            value = obj.OptionsStore.Options.IPPreprocess;
        end

        function value = get.LPMaxIter(obj)
            value = obj.OptionsStore.Options.LPMaxIter;
        end

        function value = get.LPMaxIterations(obj)
            value = obj.OptionsStore.Options.LPMaxIter;
        end

        function value = get.LPPreprocess(obj)
            value = obj.OptionsStore.Options.LPPreprocess;
        end

        function value = get.MaxNodes(obj)
            value = obj.OptionsStore.Options.MaxNodes;
        end

        function value = get.MaxNumFeasPoints(obj)
            value = obj.OptionsStore.Options.MaxNumFeasPoints;
        end

        function value = get.MaxFeasiblePoints(obj)
            value = obj.OptionsStore.Options.MaxNumFeasPoints;
        end

        function value = get.MaxTime(obj)
            value = obj.OptionsStore.Options.MaxTime;
        end

        function value = get.NodeSelection(obj)
            value = obj.OptionsStore.Options.NodeSelection;
        end

        function value = get.ObjectiveCutOff(obj)
            value = obj.OptionsStore.Options.ObjectiveCutOff;
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

        function value = get.RootLPAlgorithm(obj)
            value = obj.OptionsStore.Options.RootLPAlgorithm;
        end

        function value = get.RootLPMaxIter(obj)
            value = obj.OptionsStore.Options.RootLPMaxIter;
        end

        function value = get.RootLPMaxIterations(obj)
            value = obj.OptionsStore.Options.RootLPMaxIter;
        end

        function value = get.TolCon(obj)
            value = obj.OptionsStore.Options.TolCon;
        end

        function value = get.ConstraintTolerance(obj)
            value = obj.OptionsStore.Options.TolCon;
        end

        function value = get.RelObjThreshold(obj)
            value = obj.OptionsStore.Options.RelObjThreshold;
        end

        function value = get.ObjectiveImprovementThreshold(obj)
            value = obj.OptionsStore.Options.RelObjThreshold;
        end

        function value = get.TolFunLP(obj)
            value = obj.OptionsStore.Options.TolFunLP;
        end

        function value = get.LPOptimalityTolerance(obj)
            value = obj.OptionsStore.Options.TolFunLP;
        end

        function value = get.TolGapAbs(obj)
            value = obj.OptionsStore.Options.TolGapAbs;
        end

        function value = get.AbsoluteGapTolerance(obj)
            value = obj.OptionsStore.Options.TolGapAbs;
        end

        function value = get.TolGapRel(obj)
            value = obj.OptionsStore.Options.TolGapRel;
        end

        function value = get.RelativeGapTolerance(obj)
            value = obj.OptionsStore.Options.TolGapRel;
        end

        function value = get.TolInteger(obj)
            value = obj.OptionsStore.Options.TolInteger;
        end

        function value = get.IntegerTolerance(obj)
            value = obj.OptionsStore.Options.TolInteger;
        end
    end

    methods (Hidden)

        function thisAlgorithm = createAlgorithm(obj)
%CREATEALGORITHM Create the algorithm from the options
%
%   THISALGORITHM = CREATEALGORITHM(OBJ) creates an instance of
%   optim.algorithm.IntlinprogBranchAndCut from OBJ. The Options property
%   of THISALGORITHM is set to OBJ.

            % Create the algorithm instance and pass the options object
            thisAlgorithm = optim.algorithm.IntlinprogBranchAndCut(obj);

        end

        function checkOptions(~, ~, ~)
%CHECKOPTIONS Perform consistency checks on the options
%
%   CHECKOPTIONS(OBJ, PROBLEM, CALLER) checks that the options are
%   consistent as a set (e.g. DiffMaxChange >= DiffMinChange for fmincon).
%   This function also checks that any options that are dependent on the
%   problem are valid.

%   For this options object, there are no options that depend on eachother
%   or the problem. So this function does not need to perform any action.

        end

        function obj = setInternalOptions(obj, InternalOptions)

%SETINTERNALOPTIONS Set internal SLBI options
%
%   OBJ = SETINTERNALOPTIONS(OBJ, INTERNALOPTIONS) sets the specified
%   internal options in OBJ

           obj.InternalOptions = InternalOptions;

        end

        function OptionsStruct = extractOptionsStructure(obj)
%EXTRACTOPTIONSSTRUCTURE Extract options structure from OptionsStore
%
%   OPTIONSSTRUCT = EXTRACTOPTIONSSTRUCTURE(OBJ) extracts a plain structure
%   containing the options from obj.OptionsStore. The solver calls
%   convertForSolver, which in turn calls this method to obtain a plain
%   options structure.

            % Call the superclass method
            OptionsStruct = extractOptionsStructure@optim.options.SingleAlgorithm(obj);

            % Append the InternalOptions field
            OptionsStruct.InternalOptions = obj.InternalOptions;

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
%   the full OptionsStore. See below for an example for Intlinprog.

% Define the option defaults for the solver
OS.Defaults.BranchingRule = 'maxpscost';
OS.Defaults.CutGeneration = 'basic';
OS.Defaults.CutGenMaxIter = 10;
OS.Defaults.Display = 'iter';
OS.Defaults.HeuristicsMaxNodes = 50;
OS.Defaults.Heuristics = 'basic';
OS.Defaults.IPPreprocess = 'basic';
OS.Defaults.LPMaxIter = 'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))';
OS.Defaults.LPPreprocess = 'basic';
OS.Defaults.MaxNodes = 1e7;
OS.Defaults.MaxNumFeasPoints = Inf;
OS.Defaults.MaxTime = 7200;
OS.Defaults.NodeSelection = 'simplebestproj';
OS.Defaults.ObjectiveCutOff = Inf;
OS.Defaults.OutputFcn = [];
OS.Defaults.PlotFcns = [];
OS.Defaults.RootLPAlgorithm = 'dual-simplex';
OS.Defaults.RootLPMaxIter = 'max(30000,10*(numberOfEqualities+numberOfInequalities+numberOfVariables))';
OS.Defaults.TolCon = 1e-4;
OS.Defaults.RelObjThreshold = 0;
OS.Defaults.TolFunLP = 1e-7;
OS.Defaults.TolGapAbs = 0;
OS.Defaults.TolGapRel = 1e-4;
OS.Defaults.TolInteger = 1e-5;

% Call the package function to generate the OptionsStore
OS = optim.options.generateSingleAlgorithmOptionsStore(OS);

end

function propInfo = genPropInfo()
% Helper function to generate constant property metadata for the Intlinprog
% options class.
import optim.internal.TypeInfo;
propInfo.AbsoluteGapTolerance = TypeInfo.positiveNumericType();
propInfo.BranchRule = TypeInfo.enumType({'maxpscost','mostfractional','maxfun','strongpscost', 'reliability'});
propInfo.ConstraintTolerance = TypeInfo.numericType();
propInfo.CutGeneration = TypeInfo.enumType({'none','basic','intermediate','advanced'});
propInfo.CutMaxIterations = TypeInfo.integerType();
propInfo.Display = TypeInfo.enumType({'off','iter','final'});
propInfo.Heuristics = TypeInfo.enumType({'basic', 'intermediate','advanced','rss','rins','round','diving','rss-diving','rins-diving','round-diving','none'});
propInfo.HeuristicsMaxNodes = TypeInfo.numericType();
propInfo.IntegerPreprocess = TypeInfo.enumType({'none','basic','advanced'});
propInfo.IntegerTolerance = TypeInfo.numericType();
propInfo.LPMaxIterations = TypeInfo.integerType();
propInfo.LPOptimalityTolerance = TypeInfo.positiveNumericType();
propInfo.MaxNodes = TypeInfo.integerType();
propInfo.MaxFeasiblePoints = TypeInfo.integerType();
propInfo.MaxTime = TypeInfo.positiveNumericType();
propInfo.NodeSelection = TypeInfo.enumType({'simplebestproj','minobj','mininfeas'});
propInfo.ObjectiveCutOff = TypeInfo.numericType();
propInfo.ObjectiveImprovementThreshold = TypeInfo.positiveNumericType();
propInfo.OutputFcn = TypeInfo.fcnOrEmptyType();
propInfo.PlotFcn = TypeInfo.fcnEnumType({'optimplotmilp'});
propInfo.RelativeGapTolerance = TypeInfo.positiveNumericType();
propInfo.RootLPAlgorithm = TypeInfo.enumType({'dual-simplex','primal-simplex'});
propInfo.RootLPMaxIterations = TypeInfo.integerType();

end
