% Test of use of boundary file
%  
%  Structure of file
%    sourcefile:    OBJ file which was used to construct boundaries
%    zo_Info:       Information from TECPLOT file which was used to 
%                   associate each face to a zone
%                   zo_Info{1} = cell array of names
%                   zo_Info{2} = number of points (in TECPLOT file) assoc.
%                   with each zone
%    allBdy3D:      cell array for each zone in 3D
%                   {j}.lb: point indices for segments of boundary
%                   {j}.X,Y,Z: X,Y,Z coordinates for each segment
%                   {j}.facelist: list of faces in the zone
%   allBdy2D:      cell array for each zone in 2D
%                   {j}.lb: point indices for segments of boundary
%                   {j}.X,Y,Z: X,Y,Z coordinates for each segment
%                   {j}.facelist: list of faces in the zone
%

addpath('util_files');
addpath('util_geom');
addpath('util_plot');

%% Read geometries, if not already done
read_polyMesh_flag = 1;
read_uv_flag       = 1;
if (read_polyMesh_flag)
    % I don't love having to keep multiple copies of the polyMesh
    % and UV files. However, it's probaly better than having endless
    % confusion
    %
    
    % Replace with wherever you actually have this
    caseDir = '../caseDirs/intha_amp_50/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp fstats]=facestats([caseDir fname]);
    [R FACE S lmax] = read_polyMesh(caseDir);
end
if (read_uv_flag)
    %% Load in boundary file
    bdydir      = '../../MeshCreation/ReadTecplot_GetRegions/'
    bdyfile     = 'bdy_All_Inthavong.mat'
    
    load([bdydir bdyfile]);

    %% This will contain the proper UV file name!!!
    %% You may need to situate in your file system
    [r face r2R]=read_uv(sourcefile);

    %% Compute boundaries
    lbAll = compAllWallBdy(fp,r,face);
end

%% Need to put "olfactory boundaries" into the expected format
% Two regions: 7 and 8
Bdy3D = allBdy3D{7};
Bdy3D.lb = [Bdy3D.lb; allBdy3D{8}.lb];
Bdy3D.X  = [Bdy3D.X allBdy3D{8}.X];
Bdy3D.Y  = [Bdy3D.Y allBdy3D{8}.Y];
Bdy3D.Z  = [Bdy3D.Z allBdy3D{8}.Z];

Bdy2D = allBdy2D{7};
Bdy2D.lb = [Bdy2D.lb; allBdy2D{8}.lb];
Bdy2D.X  = [Bdy2D.X allBdy2D{8}.X];
Bdy2D.Y  = [Bdy2D.Y allBdy2D{8}.Y];
Bdy2D.Z  = [Bdy2D.Z allBdy2D{8}.Z];

% Array of time points
tarray = 0.1:0.1:3;  nT = length(tarray);


%% Do this before any plots
ucoords = r(:,1); vcoords = r(:,2);
U_patch = ucoords(face); V_patch = vcoords(face);
U_patch = U_patch'; V_patch = V_patch';

x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';


%% Create phase map if, have not already done so
create_phase_map_flag=1;
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
    Tperiod     = 3;
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
colorbar; set(gca,'FontSize',16); caxis([0,3]);
title(sprintf('Preferred time'));
linecol = cell(size(lbAll));
for k1=1:length(lbAll)
    linecol{k1} = [0 0.9 0.8];
end
optPS=[]; 
optPS.colorstr=linecol; optPS.linewidth = 1;
hold on;
plotAllWallBdy(lbAll,optPS);
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2);



if (1)
%% Just check 1 value of t: does the 2D mapping seem to make
%% sense?
    tstep = 2.25; 
    wssfname    = sprintf('%s%g/wallShearStress',caseDir,tstep);
    [w,wallInfo] = read_wss_mag_OF_face(wssfname);
    
    % Make an indicator function for region
    wIndic = zeros(size(w));
    whichFac = [allBdy3D{7}.faceList;allBdy3D{8}.faceList];
    wIndic(whichFac) = 1;
    
    
    figure;
    patch(X_patch,Y_patch,Z_patch,w,'LineStyle','None');
    axis equal;axis tight;
    colorbar; set(gca,'FontSize',16); clim([0,.25]);
    hold on;
    line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','r','LineWidth',2);
    
    
    figure;
    patch(U_patch,V_patch,w,'LineStyle','None');
    axis equal;axis tight;
    colorbar; set(gca,'FontSize',16); clim([0,.25]);
    hold on;
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','r','LineWidth',2);
end

if (0)
% Zoom in onto region.

    % LeftOlfac region
    leftLim = [0.3,0.42,0.25,0.33];
    %leftLim = [0.33,0.48,0.18,0.28];
    axis(leftLim)
    
    rightLim = [0.34,0.47,0.65,0.75];
    %rightLim = [0.33,0.48,0.74,0.84];
    axis(rightLim)
end

