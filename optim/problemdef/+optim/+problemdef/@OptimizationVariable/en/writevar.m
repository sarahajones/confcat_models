%WRITEVAR Write OptimizationVariable array to file in the standard MATLAB
%N-D grid form.
%
%   WRITEVAR(X) writes the symbolic representation of OptimizationVariable
%   array X into a text file in the standard MATLAB N-D grid form. The file
%   name is the NAME property of X, appended with '.txt'. writevar
%   overwrites any existing file.
%
%   WRITEVAR(X, FILENAME) writes the symbolic representation of an
%   OptimizationVariable array into a text file with the name given by
%   FILENAME.
%
% Examples:
%
%      % Complete use of string indexing in display:
%      x = optimvar('x', {'a','b'}, {'c','d','e'}, {'f','g','h'});
%      writevar(x, 'test.txt')
%
%      Overwrites (or creates if nonexistent) the local text file
%      'text.txt' with the following text:
% 
%      x(:,:,'f') =
% 
%         [ x('a', 'c', 'f')    x('a', 'd', 'f')    x('a', 'e', 'f') ]
%         [ x('b', 'c', 'f')    x('b', 'd', 'f')    x('b', 'e', 'f') ]
% 
%      x(:,:,'g') =
% 
%         [ x('a', 'c', 'g')    x('a', 'd', 'g')    x('a', 'e', 'g') ]
%         [ x('b', 'c', 'g')    x('b', 'd', 'g')    x('b', 'e', 'g') ]
% 
%      x(:,:,'h') =
% 
%         [ x('a', 'c', 'h')    x('a', 'd', 'h')    x('a', 'e', 'h') ]
%         [ x('b', 'c', 'h')    x('b', 'd', 'h')    x('b', 'e', 'h') ]
%
%
%      % Display an OptimizationVariable with a mix of string and numeric
%      % indicies:
%      y = optimvar('y', {'a','b'}, 3, {'f','g','h'})
%      writevar(y, 'test.txt')
%
%      Overwrites (or creates if nonexistent) the local text file
%      'text.txt' with the following text:
% 
%      y(:,:,'f') =
% 
%         [ y('a', 1, 'f')    y('a', 2, 'f')    y('a', 3, 'f') ]
%         [ y('b', 1, 'f')    y('b', 2, 'f')    y('b', 3, 'f') ]
% 
%      y(:,:,'g') =
% 
%         [ y('a', 1, 'g')    y('a', 2, 'g')    y('a', 3, 'g') ]
%         [ y('b', 1, 'g')    y('b', 2, 'g')    y('b', 3, 'g') ]
% 
%      y(:,:,'h') =
% 
%         [ y('a', 1, 'h')    y('a', 2, 'h')    y('a', 3, 'h') ]
%         [ y('b', 1, 'h')    y('b', 2, 'h')    y('b', 3, 'h') ]
%
%
%   See also OPTIMVAR, SHOWVAR, SHOWBOUNDS, WRITEBOUNDS

 
%   Copyright 2017-2018 The MathWorks, Inc.

