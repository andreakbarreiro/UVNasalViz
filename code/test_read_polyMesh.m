
% test_facestats.m
addpath('util_files');
addpath('util_geom');

whichCase= 2;
if (whichCase == 2)
    caseDir = '../caseDirs/Elad_Test_Case/';
elseif (whichCase == 1)
    % Old mesh: only 3 sided faces
    caseDir =  '/projects/abarreiro/nasoCFD/nasoCFD/intha/old_mesh/intha_normal/'

else
    % New mesh w/ boundary layers
    % Some faces have 4 sides
    caseDir =  '/projects/abarreiro/nasoCFD/nasoCFD/test_sims/intha_0.5mil_tetra/'
end

% Boundary file name
fname   = 'constant/polymesh/boundary';

%[fp fstats]=facestats([caseDir fname]);
[R FACE S lmax] = read_polyMesh(caseDir);

% We want to preserve the ability to call 
% a read function SEPARATELY from nevigating a case directory
% 

timestep    = 1.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);


%Try to plot
figure;
X=R(:,1); Y=R(:,2); Z=R(:,3);
X_patch=X(FACE)';Y_patch=Y(FACE)';Z_patch=Z(FACE)';
patch(X_patch,Y_patch,Z_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;
set(gca,'FontSize',16);


if (0)
% Try read_uv

uvfile = '../GUI/uv2.obj';
[r face r2R]=read_uv(uvfile);

% w should be defined ON THE SAME FACES
% i.e. there is a 1-1 mapping between faces in the 3D model
% and in the UV-mapping.
figure;
u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';
patch(U_patch,V_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;
set(gca,'FontSize',16);

% Read in multiple time steps
ts_array = [0.5 1 1.5 2];  nT=length(ts_array);
allWSS = cell(nT,1);
for j=1:nT
    wssfname    = sprintf('%s%g/wallShearStress',...
        caseDir,ts_array(j));
    allWSS{j} = read_wss_mag_OF_face(wssfname);
end

figure;
maxW   = 2e-4;
for j=1:nT
    w = allWSS{j};
    patch(U_patch,V_patch,w','LineStyle','None');
    axis equal;axis tight;
    colorbar;
    set(gca,'FontSize',16);
    title(sprintf('t = %g',ts_array(j)));
    caxis([0,maxW]);
    pause(1);
end

end