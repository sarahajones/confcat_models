function stop = optimplotmilp(~,optimValues,state,varargin)
%OPTIMPLOTMILP plot function for the intlinprog solver.
%   This function plots:
%   1) The objective value of the relaxed LP.
%   2) The lower bound.
%   3) The upper bound.
%   These line plots are labeled for different phases of the solver. 
%   The 'Root' legend is used to indicate both 'heuristics' and 'cutgen'
%   phases (at the root node) whereas 'Branching' legend is used for the
%   branch and bound phase.
%
%   See also intlinprog

%   Copyright 2014 The MathWorks, Inc.

% No reason to stop the solver so set the output to false.
stop = false;
% We can use either time or nodes (default) on x axis.
if nargin > 3 && ischar(varargin{1})
  useNode = ~isempty(strfind(varargin{1},'node'));
else
  useNode = true;
end
switch (state)
  case 'init'
    % Number of solution is stored in the figure appdata. Reset the appdata
    % if it exist when state = 'init'
    setappdata(gcf,'numSol',0)  
    % Legend depends on number of plots. It is cached and used when
    % plotting instead of recreating it.
    setappdata(gcf,'intlinprogLegend',{})
    hold on
    grid on
    
    if useNode
      xlabelText = 'Number of nodes';
    else
      xlabelText = 'Time elapsed';
    end
    ylabelText = 'Objective value';
    % Show x and y labels
    ylabel(ylabelText,'interp','none')
    xlabel(xlabelText,'interp','none');
    
  case {'iter','done'}
    % Information need to update the plot
    if useNode
      xdata = optimValues.numnodes;
    else  
      xdata = optimValues.time;
    end
    
    intlinprogLegend = getappdata(gcf,'intlinprogLegend');
    if optimValues.numfeaspoints > 0
      titleStr = sprintf('Best objective: %g, Relative gap: %g.', ...
        optimValues.fval,optimValues.relativegap);
    else
      titleStr = 'No feasible point yet.';
    end
    % Plot fval and lowerbound (also update the legend order if needed)
    intlinprogLegend = plotfval(optimValues,xdata,intlinprogLegend);
    intlinprogLegend = plotlowerbound(optimValues,xdata,intlinprogLegend);
    
    % Save the legend order for next call.
    setappdata(gcf,'intlinprogLegend',intlinprogLegend)    
    % Draw legend
    legend(intlinprogLegend,'Location','Best')
    title( titleStr,'interp','none');
    xlim([0,xdata+2])
end


function legendtext = plotfval(optimValues,xdata,legendtext)

% Define tags to be used in legends and line updates
tags = {'Root LB','Heuristics UB','Branch/bound UB','New Solution'};
rootlpTag = 1; heurTag = 2; bnbTag = 3; newSolTag = 4;
fval = optimValues.fval;
% Plot function value (fval)
if ~isempty(fval)
  if strcmp(optimValues.phase,'rootlp')
    % Draw marker for rootlp fval (only one call for rootlp so just a
    % point to plot)
    plot(xdata,fval,'sr','MarkerSize',10,'MarkerFaceColor','w','Tag',tags{rootlpTag});
  end

  % fval during heuristics phase
  if strcmp(optimValues.phase,'heuristics')
    plotfvalheu = findobj(get(gca,'Children'),'Tag',tags{heurTag});
    if isempty(plotfvalheu)
      plot(xdata,fval,'ok','Tag',tags{heurTag});
    else
      newX = [get(plotfvalheu,'Xdata') xdata];
      newY = [get(plotfvalheu,'Ydata') fval];
      set(plotfvalheu,'Xdata',newX, 'Ydata',newY,'Marker','o','LineStyle','-');
    end
  end
  
  % fval during branch and bound phase
  if strcmp(optimValues.phase,'branching')
    % Plot integer solution
    plotfval = findobj(get(gca,'Children'),'Tag',tags{bnbTag});
    if isempty(plotfval)
      plot(xdata,fval,'.-','Tag',tags{bnbTag});
    else
      newX = [get(plotfval,'Xdata') xdata];
      newY = [get(plotfval,'Ydata') fval];
      set(plotfval,'Xdata',newX, 'Ydata',newY);
    end
  end
  
  % Previous value of numSol
  numSol = getappdata(gcf,'numSol');
  % Draw markers for new integer solution
  if ~isempty(optimValues.numfeaspoints) && optimValues.numfeaspoints > numSol
    numSol = optimValues.numfeaspoints;
    setappdata(gcf,'numSol',numSol);
    plotnewfval = findobj(get(gca,'Children'),'Tag',tags{newSolTag});
    if isempty(plotnewfval)
      plot(xdata,fval,'*m','MarkerSize',12,'Tag',tags{newSolTag});
    else
      newX = [get(plotnewfval,'Xdata') xdata];
      newY = [get(plotnewfval,'Ydata') fval];
      set(plotnewfval,'Xdata',newX, 'Ydata',newY);
    end
  end
end

% Prepare the text for legend
for i = 1:length(tags)
  if ~any(strcmp(legendtext,tags{i})) && ...
      ~isempty(findobj(get(gca,'Children'),'Tag',tags{i}))
    legendtext{end+1} = tags{i};
  end
end


function legendtext = plotlowerbound(optimValues,xdata,legendtext)

% Define tags to be used in legends and line updates
tags = {'Branch/bound LB', 'Cuts LB'};
bnbTag = 1; cutsTag = 2;
lbound = optimValues.lowerbound;

if ~isempty(lbound)
  % Differentiate the plots between root node and b&b phase.
  if strcmp(optimValues.phase,'branching') % Branch and bound
    plotlbound_br = findobj(get(gca,'Children'),'Tag',tags{bnbTag});
    if isempty(plotlbound_br)
      plot(xdata,lbound,'.-k','Tag',tags{bnbTag});
    else
      newX = [get(plotlbound_br,'Xdata') xdata];
      newY = [get(plotlbound_br,'Ydata') lbound];
      set(plotlbound_br,'Xdata',newX, 'Ydata',newY);
    end
  else % Root node
    plotlbound = findobj(get(gca,'Children'),'Tag',tags{cutsTag});
    if isempty(plotlbound)
      plot(xdata,lbound,'or','Tag',tags{cutsTag});
    else
      newX = [get(plotlbound,'Xdata') xdata];
      newY = [get(plotlbound,'Ydata') lbound];
      set(plotlbound,'Xdata',newX, 'Ydata',newY,'Marker','o','LineStyle','-');
    end
  end
end

% Prepare the text for legend
for i = 1:length(tags)
  if ~any(strcmp(legendtext,tags{i})) && ...
      ~isempty(findobj(get(gca,'Children'),'Tag',tags{i}))
    legendtext{end+1} = tags{i};
  end
end

