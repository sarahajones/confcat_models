function optimabout()
%OPTIMABOUT helper that displays the About Box  

%   Copyright 2007-2011 The MathWorks, Inc.

a = ver('optim');
str = sprintf(['Optimization Toolbox %s\n',...
               'Copyright 1990-%s The MathWorks, Inc.'], ...
               a.Version,a.Date(end-3:end));
aboutTitle = getString(message('optim:optimtool:TitleAboutOptimTlbx'));
msgbox(str,aboutTitle,'modal');
