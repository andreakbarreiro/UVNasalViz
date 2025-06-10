%make plots for CRCNS poster
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
    
    %% Compute boundaries
    lbAll = compAllWallBdy(fp,r,face);
end

% Array of time points
tarray = 2:0.02:10;  nT = length(tarray);


%% Create phase map if, have not already done so
create_phase_map_flag=0;
if (create_phase_map_flag)
    % Set up storage for WSS
    nFace  = size(face,1);
    warray = zeros(nFace,nT);

    for j1=1:nT
        tstep = tarray(j1);
        wssfname    = sprintf('%s%g/wallShearStress',caseDir,tstep);
        [w,wallInfo] = read_wss_mag_OF_face(wssfname);
        warray(:,j1) = w;
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
end

%% Plot preferred phase map
% Get preferred phase for each point in space
% For SMOOTHED data
[prefPhaseSmo,maxOTS]   = calcPrefPhase(faceSmoAvgW,timeMod);

% use specialized colormap
RBmap = makeRedBluePhaseCmap();

figure;
patch(U_patch,V_patch,prefPhaseSmo,'LineStyle','None');
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

%% What IS the max stress where it occurs?
figure;
patch(U_patch,V_patch,maxOTS,'LineStyle','None');
colorbar; set(gca,'FontSize',16);
title(sprintf('WSS at preferred time'));
hold on;
plotAllWallBdy(lbAll,optPS);



% Cut off septal side so we can "see" turbinates
[visFACE,visface]=cutByPlane(R,FACE,r,r2R,[0 1 0 0.001]);
x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE(visFACE,:))';Y_patch=y(FACE(visFACE,:))'; 
Z_patch=z(FACE(visFACE,:))';

%% Make 3D pics of average WWS
% Peaks of ortho, retro flow
tsteps = [0.5 1.5];
for j1=1:length(tsteps)
    tind = find(timeMod==tsteps(j1));
    
    figure;
    patch(X_patch,Y_patch,Z_patch,...
        faceSmoAvgW(visFACE,tind)','LineStyle','None');
    axis equal;axis tight;
    title(sprintf('Time = %g',tsteps(j1)));
    colorbar;set(gca,'FontSize',16); caxis([0,0.25]);
    
    figure;
    patch(U_patch,V_patch,faceSmoAvgW(:,tind)','LineStyle','None');
    axis equal;axis tight;
    title(sprintf('Time = %g',tsteps(j1)));
    colorbar;set(gca,'FontSize',16); caxis([0,0.25]);
    
end

