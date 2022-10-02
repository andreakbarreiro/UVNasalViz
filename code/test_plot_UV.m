%  test w/ plotting UV
% test_plot_UV: test plotting WSS data in the UV-plane ONLY
%  You must have a previous call to read_polyMesh and read_uv.
%
addpath('util_files');
addpath('util_geom');

read_polyMesh_flag = 1;
read_uv_flag       = 1;
if (read_polyMesh_flag)
    caseDir = '../caseDirs/Elad_Test_Case_10s/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp fstats]=facestats([caseDir fname]);
    [R FACE S lmax] = read_polyMesh(caseDir);
end
if (read_uv_flag)
    uvfile = 'uv2.obj';
    [r face r2R]=read_uv([caseDir uvfile]);
end
timestep    = 8.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);

% w should be defined ON THE SAME FACES
% i.e. there is a 1-1 mapping between faces in the 3D model
% and in the UV-mapping.
figure;
u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';
patch(U_patch,V_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;


% Obtain boundary
lbAll = compAllWallBdy(fp,r,face);

hold on;
colorstr = {'r';'y';'g';'c'};
for j=1:length(lbAll)
    temp=lbAll{j};
    line(temp.Ub,temp.Vb,zeros(size(temp.Ub)),'Color',colorstr{j},'LineWidth',2);
end

set(gca,'FontSize',16);