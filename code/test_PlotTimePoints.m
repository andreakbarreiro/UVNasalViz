%
addpath('util_files');
addpath('util_geom');
addpath('util_plot');

%% Read geometries, if not already done
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
    uvfile = 'Inthavong_3D_Model.obj';
    %uvfile = 'uv2.obj';
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


FigPP = figure;
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

%%%  Now identify a few time points to plot

centroids = compCentroidsUV(r,face);

targP = [0.2 0.6;...
    0.6 0.6;...
    0.5 0.9;...
    0.5 0.82;...
    0.5 0.73;...
    0.5 0.68];

WtoPlot = cell(length(targP),1);

for j1=1:length(targP)
    % Find closest point
    nsq = sum([centroids(:,1)-targP(j1,1) centroids(:,2)-targP(j1,2)].^2,2);
    whichF = find(nsq==min(nsq));
    WtoPlot{j1} = faceSmoAvgW(whichF,:);
   
end

FigAllPts = figure;  legAllPts = cell(size(WtoPlot));
for j1=1:length(targP)
    figure; plot(timeMod,WtoPlot{j1},'linewidth',2);
    title(sprintf('u = %g, v = %g',targP(j1,1),targP(j1,2)));
    set(gca,'FontSize',16); xlabel('time');ylabel('WSS');
    figure(FigPP);hold on;plot(targP(j1,1),targP(j1,2),'*');
    plot(targP(j1,1),targP(j1,2),'ko','MarkerSize',12);
    
    figure(FigAllPts); hold on; 
    plot(timeMod,WtoPlot{j1},'linewidth',2);
    legAllPts{j1} = sprintf('u = %g, v = %g',targP(j1,1),targP(j1,2));
end

figure(FigAllPts); legend(legAllPts);
set(gca,'FontSize',16); xlabel('time');ylabel('WSS');