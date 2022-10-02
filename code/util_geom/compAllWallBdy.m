function [lbAll] = compAllWallBdy(fp,r,face)
% compAllWallBdy: compute all wall boundaries for future plotting
%
% INPUTS
%   fp:         boundary index array from call to "facestats"
%   r:          (u,v) coordinates 
%   face:       point indices for each triangular face
%
% OUTPUTS
%   lbAll:      cell structure with boundary info for each wall region
%   lball{j}.lb:   n x 2: for each bdy segment, point indices of endpoints
%   lball{j}.Ub:   2 x n: for each bdy segment, U coords of endpoints
%   lball{j}.Vb:   2 x n: for each bdy segment, V coords of endpoints
%
% Ub,Vb can be used to plot. lb can be used to investigate topology

% Get boundary
% # of wall faces per region
wallFaces = fp(:,2).*fp(:,3);
wallFaces(wallFaces==0)=[];

% Indices for each segment of the wall
startI = [0;cumsum(wallFaces)]+1;
startI = startI(1:end-1); endI = cumsum(wallFaces);

lbAll = cell(length(wallFaces),1);
for j=1:length(wallFaces)
    [lb,Ub,Vb,Zb]=findb([r zeros(length(r),1)],face(startI(j):endI(j),:));
    temp = []; temp.lb = lb; temp.Ub=Ub; temp.Vb = Vb;
    lbAll{j}=temp;
   % line(Ub,Vb,Zb,'Color',colorstr{j},'LineWidth',2);
end

end

