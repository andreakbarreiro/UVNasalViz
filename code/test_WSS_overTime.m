% test_WSS_overTime: test reading, viewing WSS data in the UV-plane over
% time
%
addpath('util_files');
addpath('util_geom');
addpath('util_plot');

%% Read geometries, if not already done
read_polyMesh_flag = 0;
read_uv_flag       = 0;
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

%% Compute boundaries
lbAll = compAllWallBdy(fp,r,face);

% Array of time points
tarray = 2:0.05:4;  nT = length(tarray);

% Set up storage for WSS
nFace  = size(face,1);
warray = zeros(nFace,nT);

for j1=1:nT
    tstep = tarray(j1);
    wssfname    = sprintf('%s%g/wallShearStress',caseDir,tstep);
    [w,wallInfo] = read_wss_mag_OF_face(wssfname);
    warray(:,j1) = w;
end

figure;
u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';

if (0)
for j1=1:nT
    patch(U_patch,V_patch,warray(:,j1)','LineStyle','None');
    axis equal;axis tight;
    colorbar; set(gca,'FontSize',16); caxis([0,0.3]);
    title(sprintf('T = %g',tarray(j1)));
    pause(0.2);
end
end

% Average over multiple periods
Tperiod     = 2;
[avgW,timeMod]      = calcAvgOverCycles(warray,tarray,Tperiod);


%% Now try smoothing in space
smoAvgW = mapQuantFaceToPoint(avgW,face,wallInfo);
% Map back to faces
faceSmoAvgW = mapQuantPointToFace(smoAvgW,face,wallInfo);


% Get preferred phase for each point in space
[prefPhase,maxOT]   = calcPrefPhase(avgW,timeMod);

% use specialized colormap
RBmap = makeRedBluePhaseCmap();

figure;
patch(U_patch,V_patch,prefPhase,'LineStyle','None');
axis equal;axis tight;colormap(RBmap);
colorbar; set(gca,'FontSize',16); caxis([0,2]);
title(sprintf('Preferred time'));

linecol = cell(size(lbAll));
for k1=1:length(lbAll)
    linecol{k1} = [0 0.9 0.8];
end
optPS=[]; 
optPS.colorstr=linecol; optPS.linewidth = 1;
hold on;
plotAllWallBdy(lbAll,optPS);


% Get preferred phase for each point in space
% For SMOOTHED data
[prefPhaseSmo,maxOTS]   = calcPrefPhase(faceSmoAvgW,timeMod);

figure;
patch(U_patch,V_patch,prefPhaseSmo,'LineStyle','None');
axis equal;axis tight;colormap(RBmap);
colorbar; set(gca,'FontSize',16); caxis([0,2]);
title(sprintf('Preferred time'));

hold on;
plotAllWallBdy(lbAll,optPS);



if (0) 
    % If you want to visualize the three possibilities
j1 = 10;
figure; subplot(1,3,1);
patch(U_patch,V_patch,avgW(:,j1)','LineStyle','None');
axis equal;axis tight;
colorbar; set(gca,'FontSize',16); caxis([0,0.3]);
title(sprintf('T = %g',tarray(j1)));

% We now use "face" in evaluating smoAvgW, because we actually need the
% values at each point on the face.
subplot(1,3,2);
temp = smoAvgW(:,j1);
patch(U_patch,V_patch,temp(face)','LineStyle','None');
axis equal;axis tight;
colorbar; set(gca,'FontSize',16); caxis([0,0.3]);

subplot(1,3,3);
patch(U_patch,V_patch,faceSmoAvgW(:,j1)','LineStyle','None');
axis equal;axis tight;
colorbar; set(gca,'FontSize',16); caxis([0,0.3]);
end