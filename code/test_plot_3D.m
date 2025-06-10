%  test w/ plotting UV
% test_plot_UV: test plotting WSS data in the UV-plane ONLY
%  You must have a previous call to read_polyMesh and read_uv.
%
addpath('util_files');
addpath('util_geom');

read_polyMesh_flag = 0;
read_uv_flag       = 0;
if (read_polyMesh_flag)
    caseDir = '../caseDirs/Inthavong_OrthoOnly_0p5Hz/';
    %caseDir = '../caseDirs/Elad_Test_Case_10s/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp fstats]=facestats([caseDir fname]);
    [R FACE S lmax] = read_polyMesh(caseDir);
end
if (read_uv_flag)
    uvfile = 'uv2.obj';
    [r face r2R]=read_uv([caseDir uvfile]);
end
timestep    = 0.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);

figure;
x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';
patch(X_patch,Y_patch,Z_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;

if (0)
% Does "getplane" do what I think it does??
[visFACE,visface]=cutByPlane(R,FACE,r,r2R,[0 1 0 0.001]);

% ONLY plot visible faces
figure;
X_patch=x(FACE(visFACE,:))';Y_patch=y(FACE(visFACE,:))'; 
Z_patch=z(FACE(visFACE,:))';
patch(X_patch,Y_patch,Z_patch,w(visFACE)','LineStyle','None');
axis equal;axis tight;
colorbar;
end

if (0)
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

end