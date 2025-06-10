function cmap = makeTwoColorCmap(c1,c2,myopt)
% Make c1 <-> c2 colormap
%
% myopt: include shading? [0/1] default = 0

shadeType = 0;
if (isfield(myopt,'shadeType')); shadeType = myopt.shadeType; end;

% Linear interpolate between c1 and white; white and c2
dh = 0.025;

xval = 0:dh:1; xval = xval';

%rmap = [ones(size(xval)) xval xval];

c1map = (1-xval)*c1 + xval*[1 1 1];

%bmap = [xval xval ones(size(xval))];
c2map = (1-xval)*c2 + xval*[1 1 1];
c2map = flipud(c2map);

cmap = [c1map;c2map(2:end,:)];

if (shadeType == 1)
    nxval = xval.^4;
    nxval = [nxval; flipud(nxval(2:end))];
    cmap = cmap .* repmat(nxval,1,3);
end