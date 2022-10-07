function RBmap = makeRedBluePhaseCmap()
% Make red <-> blue colormap
% suitable for ortho/retro PHASE 
dh = 0.02;

xval = 0:dh:1; xval = xval';
nL   = length(xval);

if (0)
% Cut down on possible saturation?
dsat = 0.2;
xval = dsat + xval*(1-dsat);
end

% Make this double sided
% We get "white areas" towards the end of the 
%   breath cycle
xval = [flipud(xval);xval(2:end)];

rmap = [ones(size(xval)) xval xval];
bmap = [xval xval ones(size(xval))];

% Put blue before red
rmap = flipud(rmap);

% Introduce gray-scaling to get early vs. late cycle?
rmap(nL+1:end,:) = rmap(nL+1:end,:)*0.8;
bmap(nL+1:end,:) = bmap(nL+1:end,:)*0.8;

RBmap = [bmap;rmap(2:end,:)];