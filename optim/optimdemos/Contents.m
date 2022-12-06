Optimization Toolbox 
%
% Demonstrations.
%   portfoptimdemo - Quadratic programming to solve portfolio optimization problem.
%   circustent   - Quadratic programming to find shape of a circus tent.
%   optdeblur    - Image deblurring using bounded linear least-squares.
%   symbolic_optim_demo - Use symbolic toolbox functions to compute
%                         gradients and Hessians.
%   tutdemo      - Tutorial walk-through.
%   optimparfor  - Minimize an expensive optimization problem im parallel (requires
%                  Parallel Computing Toolbox and Global Optimization Toolbox).
%   goaldemo     - Goal attainment.
%   dfildemo     - Finite-precision filter design (requires Signal Processing
%                  Toolbox).
%   datdemo      - Fitting a curve to data.
%   officeassign - Binary integer programming to solve the office assignment
%                  problem.
%   bandem       - Banana function minimization demonstration.
%   airpollution - Use semi-infinite programming to analyze the effect of
%                  uncertainty.
%   SudokuExample - Sudoku problems using integer linear programming.
%   TravellingSalesmanExample - Travelling salesman problem using integer
%                  linear programming.
%   FactoryExample - Optimal production and distribution levels using integer
%                  linear programming.
%   OptimalDispatchExample - Optimal scheduling of power generators using
%                  integer linear programming.
%   PortfolioMIQPExample - Mixed-integer quadratic programming portfolio
%                  optimization.
%   LongTermInvestmentsExample - Optimize investments over a time period
%                  using linear programming.
% 
% Examples from User's Guide
%   objfun       - nonlinear objective
%   confun       - nonlinear constraints
%   objfungrad   - nonlinear objective with gradient
%   confungrad   - nonlinear constraints with gradients
%   confuneq     - nonlinear equality constraints
%   optsim.mdl   - Simulink model of nonlinear plant process
%   optsiminit   - parameter initialization for optsim.mdl
%   bowlpeakfun  - objective function for parameter passing
%   nestedbowlpeak - nested objective function for parameter passing
%   nlsf1         - nonlinear equations objective with Jacobian
%   nlsf1a        - nonlinear equations objective 
%   nlsdat1       - MAT-file of Jacobian sparsity pattern (see nlsf1a)
%   brownfgh      - nonlinear minimization objective with gradient and Hessian
%   brownfg       - nonlinear minimization objective with gradient 
%   brownhstr     - MAT-file of Hessian sparsity pattern (see brownfg)
%   browneq       - MAT-file of Aeq and beq sparse linear equality constraints
%   runfleq1      - demonstrates 'HessMult' option for FMINCON with equalities
%   brownvv       - nonlinear minimization with dense structured Hessian
%   hmfleq1       - Hessian matrix product for brownvv objective
%   fleq1         - MAT-file of V, Aeq, and beq for brownvv and hmfleq1 
%   qpbox1        - MAT-file of quadratic objective Hessian sparse matrix
%   runqpbox4     - demonstrates 'HessMult' option for QUADPROG with bounds
%   runqpbox4prec - demonstrates 'HessMult' and TolPCG options for QUADPROG
%   qpbox4        - MAT-file of quadratic programming problem matrices
%   particle      - MAT-file of linear least squares C and d sparse matrices
%   sc50b         - MAT-file of linear programming example
%   densecolumns  - MAT-file of linear programming example
%

% Internally Used Utility Routines
%
%   Demonstration utility routines
%   elimone           - eliminates a variable (used by dfildemo)
%   filtobj           - frequency response norm (used by dfildemo)
%   filtcon           - frequency response roots (used by dfildemo)
%   fitvector         - value of fitting function (used by datdemo)
%   tentdata          - MAT-file of data for circustent demo
%   optdeblur         - MAT-file of data for optdeblur demo
%   plotdatapoints    - helper plotting function (used by datdemo)
%   printofficeassign - helper plotting function (used by officeassign demo)
%   filtfun           - returns frequency response and roots (used by dfildemo)
%   filtfun2          - returns frequency response norm and roots (used by dfildemo) 
%   plotPortfDemoStandardModel.m - generates plots (used by portfoptimdemo)
%   plotPortfDemoGroupModel.m - generates plots (used by portfoptimdemo)
%   detectSubtours.m  - returns subtours (used by TravellingSalesmanExample)
%   updateSalesmanPlot.m - helper plotting function (used by TravellingSalesmanExample)
%   drawSudoku.m      - helper plotting function (used by SudokuExample)
%   sudokuEngine.m    - solves Sudoku puzzles (used by SudokuExample)
%   dispatchPrice     - MAT-file of data for OptimalDispatchExample
%   plotInvestments   - generates plots (used by LongTermInvestmentsExample)
%   eil33.2.mps       - data file for mpsread (used by function reference page)
%

%   Copyright 1990-2018 The MathWorks, Inc.
