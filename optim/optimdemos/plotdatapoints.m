function h = plotdatapoints(t,y)
%PLOTDATAPOINTS Helper function for DATDEMO

%   Copyright 1990-2008 The MathWorks, Inc.

h = plot(t,y,'b-');
axis([0 2 -0.5 6])
hold on
plot(t,y,'ro')    
title('Data points and fitted curve')
hold off
