function err = optimguiGeneratemfile(hashProb,hashOpt)
%optimguiGeneratemfile generates a MATLAB file from OPTIMTOOL.
%   hashProb and hashOpt are Java hash tables containing information about
%   the problem and options model. hashProb and hashOpt contain only information 
%   that user has changed since last time the data (Java model) from the GUI 
%   was passed to MATLAB workspace. (E.g. at the time of exporting, running,
%   generating code.) This function will update the MATLAB workspace and will 
%   call generateMfile.m to generate a MATLAB file.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2015 The MathWorks, Inc.

err = '';
% Get modified fields from the GUI
[probStruct,optStruct,errProb,errOpt] = readOptimHashTable(hashProb, hashOpt);
if ~isempty(errProb)
    err = sprintf('%s\n',getString(message('optim:optimtool:RuntimeErrorHeaderProblem',errProb)));
    return;
elseif ~isempty(errOpt)
    err = sprintf('%s\n',getString(message('optim:optimtool:RuntimeErrorHeaderOptions',errOpt)));
    return;
end
optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
if ~isempty(optimtoolGui)
    currentoptFields = fieldnames(optStruct);
    for i = 1:length(currentoptFields)
        if ~javaMethodEDT('isOptionEnabled',optimtoolGui,currentoptFields(i))
            % Remove the disabled option from the optStruct (not needed to
            % generate MATLAB file)
            optStruct = rmfield(optStruct,currentoptFields{i});
        end
    end
end

% Generate MATLAB file with for modified problem and options structure
generateMfile(probStruct,optStruct);

%----------------------------------------------------------------------
function msg = generateMfile(probStruct,optStruct)
%generateMfile helper function for optimguiGeneratemfile. 
%
%   This function takes a problem structure 'probStruct' and an options
%   structure 'optStruct' and generates a MATLAB file. The generated code sets
%   the non-default options and also calls the appropriate solver.
%

msg  = '';
% headerCode is the code for the generated function signature
headerCode = 'function [';
% optionsCode is the code for options settings
optionsCode = '';
% solverCode is the code for calling the appropriate solver
solverCode = '[';

% Get file name to use, remember the directory name
filespec = '*.m';
[optimMfileName,pn] = uiputfile(filespec, ...
    getString(message('optim:optimtool:PutFileDlgTitleGenFile')),'untitled.m');
if isequal(optimMfileName,0) || isequal(pn,0)
    return
end
if ~ismember('.',optimMfileName)
    optimMfileName = [optimMfileName '.m'];
end
optimMfileName = sprintf('%s%s',pn,optimMfileName);


% Get MATLAB file name with .m suffix, and get corresponding function name
if length(optimMfileName)<2 || ~isequal(optimMfileName(end-1:end),'.m')
    optimMfileName = sprintf('%s.m',optimMfileName);
end
[dirname,fcnname] = fileparts(optimMfileName);

solver = probStruct.solver; 
% Problem structure will tell us about the number of input arguments and
% also about the order in which the solver accept them. The input 'probStruct'
% contains fields for all the solvers. The output 'probStruct' will only 
% contain fields relevant to 'solver'.
probStruct = createProblemStruct(solver,[],probStruct);
% We do not want the 'solver' field in the problem structure
probStruct = rmfield(probStruct,'solver');
% Check if random states are used in MATLAB file generation
if isfield(probStruct,'rngstate')
    probStruct = rmfield(probStruct,'rngstate');
end

problemFields = fieldnames(probStruct);
problemNumFields = length(problemFields);

% Result structure will tell us about the number of output arguments and
% also about the order in which it comes from the solver.
resultStruct = createResultsStruct(solver);
resultsFields = fieldnames(resultStruct);
resultsNumFields = length(resultsFields);
% Append solver's output list to the  signature of the generated 
% function and also to the calling syntax for the solver.
for i = 1:resultsNumFields
   headerCode = [headerCode sprintf('%s,',resultsFields{i})]; 
   solverCode = [solverCode sprintf('%s,',resultsFields{i})]; 
end
% Replace the last char ',' (comma) by ']' in the headerCode/solverCode 
headerCode(end) = ']'; solverCode(end) = ']';
% Add '=' character to the end
headerCode = [headerCode ' = ']; solverCode = [solverCode ' = '];
% Add the file name to headerCode
headerCode = [headerCode sprintf('%s(',fcnname)];
% Add the 'solver' name to solverCode
solverCode = [solverCode sprintf('...\n') sprintf('%s(',solver)];

% For each field of the problem structure
for i = 1:problemNumFields
    probField = problemFields{i};
    probValue = probStruct.(probField);
    tempcode = getStringValue(probValue,probField,true);
    solverCode = [solverCode tempcode sprintf(',')];
        % If value is numeric and non-empty then this code also goes into
        % the headerCode.
        if ~isempty(probValue)  && isnumeric(probValue)
            headerCode = [headerCode tempcode sprintf(',')];  
        end
end

% After all the problem data is written we want to add options as the last
% argument to solverCode
solverCode = [solverCode sprintf('options);')];
% Get fieldnames of the options structrue
[optStruct, optionsFcn] = createOptionsStruct(solver,optStruct);
% optStruct should have only the valid fields for optimoptions.
[~, optStruct] = optimset2optimoptions(solver,optStruct);
optionsFields = fieldnames(optStruct);
optionsNumFields = length(optionsFields);

% Start with default options
optionsCode = [optionsCode getString(message('optim:optimtool:CommentDefOpts'))];
if strcmp(optionsFcn,'optimoptions')
    optionsCode = [optionsCode sprintf('options = %s(''%s'');\n',optionsFcn,solver)];
    % Create an options object from the structure
    optionsClassName = sprintf('%s%s', upper(solver(1)), solver(2:end));
    if any(strcmp(solver, {'ga', 'gamultiobj', 'patternsearch', 'simulannealbnd'}))
        optionsClassName = [optionsClassName, 'Options'];
    end
    thisOptions = optim.options.(optionsClassName);    
else
    optionsCode = [optionsCode sprintf('options = %s;\n',optionsFcn)];
    thisOptions = [];
end

% For each property
optionsCode = [optionsCode getString(message('optim:optimtool:CommentModifyOpts'))];
for i = 1:optionsNumFields
    optField = optionsFields{i};
    optValue = optStruct.(optField);
    
    % Don't generate code for defaults or when 
    
    % 1. Hessian is set to 'on' or 'user-supplied'.    
    % This is because Hessian maps to HessianApproximation and we cannot
    % set HessianApproximation to 'on' or 'user-defined'. In this case the
    % value of HessianApproximation is ignored (setting Hessian to
    % 'on'/'user-defined' means that a HessianFcn or HessianMultiplyFcn is
    % going to be specified), so it is fine not to generate code for this
    % case.
        
    % 2. FitnessLimit or StallTimeLimit is set for gamultiobj 
    % optimtool has left options widgets for FitnessLimit and
    % StallTimeLimit for gamultiobj. However, these options are not
    % supported by gamultiobj. As we're planning to deprecate optimtool,
    % we'll leave the options in the tool and just not generate code here.
    if  ~isempty(optValue) && ...
            ~(strcmp(optField, 'Hessian') && any(strcmp(optValue, {'on', 'user-supplied'}))) && ...
            ~(strcmpi(optField, 'fitnesslimit') && strcmpi(solver, 'gamultiobj')) && ...
            ~(strcmpi(optField, 'stalltimelimit') && strcmpi(solver, 'gamultiobj')) 
        
        if isempty(thisOptions)
            optFieldNewName = optField;
            optFieldNewValue = optValue;
        else
            % Set the old property in the options object. This will allow us to
            % get the correct new value using the mapping code in the options
            % object get methods.
            thisOptions.(optField) = optValue;
            % Get the new name
            optFieldNewName = optim.options.OptionAliasStore.getNameFromAlias(optField);
            optFieldNewName = optFieldNewName{1};
            
            % If the new name is TolFun and the solver is one of the global
            % solvers, then the new name should be FunctionTolerance
            if strcmp(optField, 'TolFun') && ...
                    any(strcmp(solver, {'ga', 'gamultiobj', 'patternsearch', 'simulannealbnd'}))
                optFieldNewName = 'FunctionTolerance';
            end
            
            % Get the new value. If the new value doesn't exist in the solver
            % options object, use the existing name and value
            try
                optFieldNewValue = thisOptions.(optFieldNewName);
            catch
                optFieldNewName = optField;
                optFieldNewValue = optValue;
            end
        end
        tempcode = getStringValue(optFieldNewValue,optFieldNewName,false);
        optionsCode = [optionsCode sprintf('options = %s(options,''%s'', %s);\n', ...
            optionsFcn,optFieldNewName,tempcode)];
        
        % If TolFun has been set in optimtool, we need to set
        % OptimalityTolerance and FunctionTolerance for those algorithms
        % that have both of these options.
        if ~isempty(thisOptions) && strcmp(optFieldNewName, 'OptimalityTolerance')
            try %#ok
                thisOptions.FunctionTolerance;
                optionsCode = [optionsCode sprintf(...
                    'options = optimoptions(options,''FunctionTolerance'', %s);\n', ...
                    tempcode)];
            end
        end
      
        % If value is numeric and non-empty then this code also goes into
        % the headerCode.
        if ~isempty(optFieldNewValue) && isnumeric(optFieldNewValue)
            headerCode = [headerCode tempcode sprintf(',')];    
        end
    end
end

% Replace the last char ',' (comma) by the ')' (closing bracket) in the headerCode
if headerCode(end) == ','
    headerCode(end) = ')';
else
    headerCode(end) = '';
end
% Add some generic comments in the file
headerCode = [headerCode getString(message('optim:optimtool:CommentMfileHeader'))];
%headerCode = [headerCode sprintf('%s\n','% Optimization Toolbox is required to run this MATLAB file.')];

% Concatenate all the generated code into one
code = [sprintf('%s',headerCode),sprintf('%s',optionsCode),sprintf('%s',solverCode)];
% Open and write to the file
[fid,fopenMsg] = fopen(optimMfileName,'w');
if fid==-1
    msg = getString(message('optim:optimtool:ErrorGenMfile',optimMfileName,fopenMsg));
    errordlg(msg,getString(message('optim:optimtool:ErrorDlgTitleCantSave')),'modal');
    return
end
% Close the file
fprintf(fid,'%s\n',code);
st = fclose(fid);
if st ~= 0
    msg = getString(message('optim:optimtool:ErrorClosingFile',fcnname));
    return;
end
% Open the MATLAB file just created
edit(optimMfileName)

function rhsValue = getStringValue(value,FieldName,problem_type)
% Function to wrap around value2RHS so that matrices are not
% converted to strings. 
% If problem_type is 'true' then string for solver's input  arguments are generated
% otherwise for options using the 'FieldName'

if ~isempty(value) && isnumeric(value)
    if problem_type
        rhsValue = sprintf('%s',FieldName); % This is for problem values that are non-empty
    else
        rhsValue = sprintf('%s_Data',FieldName); % This is for options values that are matrices/vectors
    end
elseif islogical(value)
    % value2RHS does handle the logical case. However it maps it to the
    % string '1' or '0', which is what optimtool requires. We require the
    % strings 'true' or 'false' for code generation.
    if value
        rhsValue = 'true';
    else
        rhsValue = 'false';
    end    
else
    rhsValue = value2RHS(value);
end

