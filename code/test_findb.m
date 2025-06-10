
% test_facestats.m
addpath('util_files');
addpath('util_geom');

caseDir = '../Elad_Test_Case/';

% Boundary file name
fname   = 'constant/polymesh/boundary';

[fp fstats]=facestats([caseDir fname]);
[R FACE S lmax] = read_polyMesh(caseDir);

% # of wall faces per region
wallFaces = fp(:,2).*fp(:,3);
wallFaces(wallFaces==0)=[];

startI = [0;cumsum(wallFaces)]+1;
startI = startI(1:end-1); endI = cumsum(wallFaces);

figure;hold on;
for j=1:4
    [lb,Xb,Yb,Zb]=findb(R,FACE(startI(j):endI(j),:));

    line(Xb,Yb,Zb,'Color','r','LineWidth',2);
end
% % We want to preserve the ability to call 
% % a read function SEPARATELY from nevigating a case directory
% % 
% 
% timestep    = 1.5;
% wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);
% 
% w = read_wss_mag_OF_face(wssfname);
