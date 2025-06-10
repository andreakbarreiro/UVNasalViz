function RBmap = makeHSVPhaseCmap()
% Make HSV colormap
% suitable for ortho/retro PHASE 

temp = colormap(hsv);

% Output has 201 colors ranging from red->orange ->...

% This is an HSV map with the property that 
% Blue sits at 0.25
% Red sits at 0.75
RBmap = [temp(101:end,:);temp(1:34,:)];
