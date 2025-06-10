% make plots for CRC Talk
% April 2023
%


addpath('util_files');
addpath('util_geom');
addpath('util_plot');

%% Read geometries, if not already done
read_polyMesh_flag = 0;
read_uv_flag       = 0;
if (read_polyMesh_flag)
    % I don't love having to keep multiple copies of the polyMesh
    % and UV files. However, it's probaly better than having endless
    % confusion
    %
    % Inthavong directory    
    %caseDir = '../caseDirs/intha_amp_50/';
    caseDir =  '/projects/abarreiro/nasoCFD/nasoCFD/intha/old_mesh/intha_normal/'

% /intha/old_mesh
%   /intha_normal
%   /intha_lam_freq_0p17
%   /intha_lam_freq_0p64
%   /intha_totpress_normal
% /intha/new_mesh

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp fstats]=facestats([caseDir fname]);
    [R FACE S lmax] = read_polyMesh(caseDir);
end
if (read_uv_flag)
    %% Load in boundary file
    %bdydir      = '../../MeshCreation/ReadTecplot_GetRegions/'

    % On M3
    bdydir = '../caseDirs/Inthavong_Test_Case/'
    bdyfile     = 'bdy_All_Inthavong.mat'
    %bdyfile     = 'bdy_Olfac_Inthavong.mat';
    load([bdydir bdyfile]);

    %% This will contain the proper UV file name!!!
    [r face r2R]=read_uv(sourcefile);
    if (0)
        if (0)
            % Issue: UV coordinates are reversed, somehow, from uv2.obj
            % 
            % We appear to have the same number of faces, however. WHY???
            %
            uvfile = 'uv.obj';
            [r face r2R]=read_uv([caseDir uvfile]);
        end
        objdir  = '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/';
        objfile = 'uv2.obj';
    
        [r face r2R]=read_uv([objdir objfile]);
    end

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
tarray = 0.02:0.02:3;  nT = length(tarray);


%% Do this before any plots
ucoords = r(:,1); vcoords = r(:,2);

U_patch = ucoords(face); V_patch = vcoords(face);
U_patch = U_patch'; V_patch = V_patch';

x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';


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

if (0)
%% Now, I want to find a specific point to plot stress vs. time
Umean=mean(U_patch); Vmean=mean(V_patch);
mydist = (Umean-0.2).^2 + (Vmean-0.2).^2;
whichind = find(mydist==min(mydist),1);

%% I am not sure what solution I want for this.
%% This is just a stopgap
timeMod(end)=Tperiod;

figure;
plot(timeMod,avgW(whichind,:),'.-','linewidth',2);
hold on;set(gca,'FontSize',16); 
%plot(timeMod,faceSmoAvgW(whichind,:),'.-');
end

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

if (1)
    % Load in Abdullah's 2D figures, Zoom and view

%% Orthonasal peak
figloc = '../MatLab_figs/'
wssOrthPeakFig = openfig([figloc 'wss075.fig']);

%% Retronasal peak
wssRetPeakFig = openfig([figloc 'wss225.fig']);

%% Preferred time
wssPrefFig = openfig([figloc 'prefTime.fig']);

figure(wssOrthPeakFig)
clim([0,0.3]);
hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','r','LineWidth',2);

figure(wssRetPeakFig)
clim([0,0.3]);
hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','r','LineWidth',2);

% LeftOlfac region
leftLim = [0.3,0.42,0.25,0.33];
%leftLim = [0.33,0.48,0.18,0.28];
axis(leftLim)

rightLim = [0.34,0.47,0.65,0.75];
%rightLim = [0.33,0.48,0.74,0.84];
axis(rightLim)

figure(wssPrefFig)
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2);

end

if (1)
%% 3D figures
figloc = '../MatLab_figs/'
wssOrthPeakFig = openfig([figloc 'wss075_3d.fig']);

figure(wssOrthPeakFig)
hold on;
line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','r','LineWidth',2);

view(55,42)
clim([0,0.3])
set(gca,'FontSize',16)

wssRetPeakFig = openfig([figloc 'wss225_3d.fig']);
figure(wssRetPeakFig);
hold on;
line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','r','LineWidth',2);

view(55,42)
clim([0,0.3])
set(gca,'FontSize',16)

end