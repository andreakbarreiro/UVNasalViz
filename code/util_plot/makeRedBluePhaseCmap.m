function RBmap = makeRedBluePhaseCmap()
% Make red <-> blue colormap
% suitable for ortho/retro PHASE 
dh = 0.05;

xval = 0:dh:1; xval = xval';
nL   = length(xval);

% For each segment we interpolate between two colors
%
% BLUE
endA1 = [1 1 1];
endA2 = [0 0 1];
bmap = (1-xval)*endA1 + xval*endA2;

% VIOLET
endB1 = endA1;
endB2 = [0.8 0 1];
vmap = (1-xval)*endB1 + xval*endB2;

% RED
endC1 = endA1;
endC2 = [1 0 0];
rmap = (1-xval)*endC1 + xval*endC2;

% ORANGE
endD1 = endA1;
endD2 = [1 0.7 0];
omap = (1-xval)*endD1 + xval*endD2;

RBmap = [bmap; flipud(vmap(2:end,:));...
    rmap; flipud(omap(2:end,:))];

if (0)
% Make this double sided
% We get "white areas" towards the end of the 
%   breath cycle
xval = [flipud(xval);xval(2:end)];

rmap = [ones(size(xval)) xval xval];
bmap = [xval xval ones(size(xval))];

% Put blue before red
rmap = flipud(rmap);

% Add a second color, to give more dimension to the cycle?
% Red ->orange
rmap(nL+1:end,2) = rmap(nL+1:end,1)*0.9;

% Blue -> purple
bmap(nL+1:end,1) = bmap(nL+1:end,3)*0.9;

% Introduce gray-scaling to get early vs. late cycle?
%rmap(nL+1:end,:) = rmap(nL+1:end,:)*0.8;
%bmap(nL+1:end,:) = bmap(nL+1:end,:)*0.8;

RBmap = [bmap;rmap(2:end,:)];

end