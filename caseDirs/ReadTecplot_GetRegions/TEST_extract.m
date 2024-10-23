%% Now we visualize and test


addpath('../../UVViz/code/util_files');
addpath('../../UVViz/code/util_geom');

% This should be the boundary file
bdyfile    = 'bdy_All_Inthavong.mat';
load(bdyfile);

% The boundary file will have the OBJ file
% Read in obj file
[r,face,r2R,R]=read_3D_from_OBJ(sourcefile);



% Where to find a WSS file
caseDir = '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/'
timestep    = 0.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);
FACE = r2R(face);
x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';

u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';

f3d = figure;
f2d = figure;

nReg = length(zo_Info{1});

for j1=1:nReg
    figure(f3d);  
    patch(X_patch,Y_patch,Z_patch,w','LineStyle','None');
    axis equal;axis tight; colorbar;
    hold on;
    X3 = allBdy3D{j1}.X; Y3 = allBdy3D{j1}.Y; Z3 = allBdy3D{j1}.Z;
    line(X3,Y3,Z3,'Color','r','LineWidth',2);
    set(gca,'FontSize',16);

    figure(f2d);  
    patch(U_patch,V_patch,w','LineStyle','None');
    axis equal;axis tight; colorbar;

    hold on;
    X2 = allBdy2D{j1}.X; Y2 = allBdy2D{j1}.Y; Z2 = allBdy2D{j1}.Z;
    line(X2,Y2,Z2,'Color','r','LineWidth',2);

    zo_Info{1}{j1}
    pause;

    figure(f3d); clf;
    figure(f2d); clf;
   
end
