function [ax] = plotAllWallBdy(lbAll,optPS)
% plot all of the wall boundaries
% 
% INPUTS
%   lbAll:      cell structure with boundary info for each wall region
%   lball{j}.lb:   n x 2: for each bdy segment, point indices of endpoints
%   lball{j}.Ub:   2 x n: for each bdy segment, U coords of endpoints
%   lball{j}.Vb:   2 x n: for each bdy segment, V coords of endpoints
%    
% INPUTS, OPTIONAL:
%   optPS.colorstr:   color specifiers for different parts of wall bdy

% DEFAULTS
colorstr = cell(length(lbAll),1);
for j1=1:length(lbAll)
    colorstr{j1}='r';
end
lw       = 2;

% Process optional arguments
if (isfield(optPS,'colorstr')); colorstr = optPS.colorstr; end 
if (isfield(optPS,'linewidth')); lw = optPS.linewidth; end 
      

for j=1:length(lbAll)
    temp=lbAll{j};
    line(temp.Ub,temp.Vb,zeros(size(temp.Ub)),'Color',colorstr{j},...
        'LineWidth',lw);
end


end

