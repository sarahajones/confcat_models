function colour = mT_pickColour(colourNum)
% Returns 3-element vector specifying a colour. Each input colourNum, provides a
% different output colour.

% Joshua Calder-Travis, j.calder.travis@gmail.com

colours = NaN(7, 3);

% Colours taken from Ito & Okabe, Colour Universal Design
% (https://jfly.uni-koeln.de/color/index.html)
colours(1, :) = [0, 0, 0];
colours(2, :) = [0.9, 0.6, 0];
colours(3, :) = [0.35, 0.7, 0.9];
colours(4, :) = [0, 0.6, 0.5];
colours(5, :) = [0, 0.45, 0.7];
colours(6, :) = [0.8, 0.4, 0];
colours(7, :) = [0.8, 0.6, 0.7];

% Pick!
colourNum = mod(colourNum, 7);
colourNum(colourNum == 0) = 7;
colour = colours(colourNum, :);