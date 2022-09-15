%  test w/ plotting UV
% test_plot_UV: test plotting WSS data in the UV-plane ONLY
%  You must have a previous call to read_polyMesh and read_uv.
%
addpath('util_files');
addpath('util_geom');

read_polyMesh_flag = 1;
read_uv_flag       = 1;
if (read_polyMesh_flag)
    caseDir = '../caseDirs/Elad_Test_Case/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    %[fp fstats]=facestats([caseDir fname]);
    [R FACE S lmax] = read_polyMesh(caseDir);
end
if (read_uv_flag)
    uvfile = 'uv2.obj';
    [r face r2R]=read_uv([caseDir uvfile]);
end
timestep    = 1.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);

% Get boundary
% # of wall faces per region
wallFaces = fp(:,2).*fp(:,3);
wallFaces(wallFaces==0)=[];

startI = [0;cumsum(wallFaces)]+1;
startI = startI(1:end-1); endI = cumsum(wallFaces);


% w should be defined ON THE SAME FACES
% i.e. there is a 1-1 mapping between faces in the 3D model
% and in the UV-mapping.
figure;
u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';
patch(U_patch,V_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;

hold on;
lbAll = cell(length(wallFaces),1);
for j=1:length(wallFaces)
    [lb,Ub,Vb,Zb]=findb([r zeros(length(r),1)],face(startI(j):endI(j),:));
    temp = []; temp.lb = lb; temp.Ub=Ub; temp.Vb = Vb;
    lbAll{j}=temp;
    %line(Ub,Vb,Zb,'Color','r','LineWidth',2);
end

set(gca,'FontSize',16);