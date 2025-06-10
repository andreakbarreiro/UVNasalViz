% Compare phase preferences
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
faceList = [allBdy2D{7}.faceList; allBdy2D{8}.faceList];

%% Do this before any plots
ucoords = r(:,1); vcoords = r(:,2);

U_patch = ucoords(face); V_patch = vcoords(face);
U_patch = U_patch'; V_patch = V_patch';

x=R(:,1); y=R(:,2);  z=R(:,3);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';


%% Do not allow creation of phase mape here.
%% Use calcPrefPhase_driver



% use specialized colormap
RBmap = makeRedBluePhaseCmap();

compare_two_periods = 1;
compare_two_options = 0;

% Two different periods
if (compare_two_periods)
    %% Load pre-computed maps
    dirlist = {'../caseDirs/intha_amp_50_T6/';...
        '../caseDirs/intha_amp_50_T1p5/';...
        };
    %% Assume name is the same in each case
    filelist{1} = [dirlist{1} 'prefPhaseSmo_retroAdjust.mat'];
    filelist{2} = [dirlist{2} 'prefPhaseSmo_retroAdjust.mat'];
elseif (compare_two_options)
    %% Load pre-computed maps
    dirlist = {'../caseDirs/intha_amp_50/';...
        '../caseDirs/intha_amp_50/';...
        };
    %% Assume name is the same in each case
    filelist{1} = [dirlist{1} 'prefPhaseSmo.mat'];
    filelist{2} = [dirlist{2} 'prefPhaseSmo_retroAdjust.mat'];
end

prefPhaseSmoCell = cell(length(dirlist),1);
Tarray    = zeros(length(dirlist),1);

for j1=1:length(dirlist)
    temp=load(filelist{j1});
    prefPhaseSmoCell{j1}=temp.prefPhaseSmo;
    Tarray(j1) = temp.Tperiod;
end

for j1=1:length(dirlist)
    % do I need to normalize?
    prefPhaseSmoCell{j1}=prefPhaseSmoCell{j1}/Tarray(j1);
end

if (0)
    % Do we normalize the "1.5" thing?
    prefPhaseSmoCell{2}=prefPhaseSmoCell{2}*1.5/(1/0.68);

end


if (compare_two_periods)
    title_str_cell = {sprintf('T = %g',Tarray(1)),...
    sprintf('T = %g',Tarray(2))};
else
    title_str_cell = {sprintf('T = %g, no adjustment',Tarray(1)),...
    sprintf('T = %g, retro-adjusted ',Tarray(2))};
end

% Each preferred phase in grayscale
figure; 
for j1=1:2
    subplot(1,2,j1)
    patch(U_patch,V_patch,prefPhaseSmoCell{j1},'LineStyle','None');
    axis equal;axis tight;colormap(RBmap);
    colorbar; set(gca,'FontSize',16); caxis([0,1]);
    %title(sprintf('T = %g',Tarray(j1)));
    title(title_str_cell{j1});
    hold on;
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2);
end


myopt = struct('shadeType',0);
hsvTrunc = makeTwoColorCmap([0 1 0],[1 0 1],myopt);

%% The difference
mydiff = prefPhaseSmoCell{2}-prefPhaseSmoCell{1};

% Do a threshold of significance
sigVal = 0.012;
mydiffThresh = (abs(mydiff)>=sigVal).*sign(mydiff);

%% Thresholded to "1" or "-1" if above significance
figure;subplot(1,2,1);
%patch(U_patch,V_patch,prefPhaseSmoCell{2}-prefPhaseSmoCell{1},'LineStyle','None');
patch(U_patch,V_patch,mydiffThresh,'LineStyle','None');
axis equal;axis tight;colormap(hsvTrunc);
colorbar; set(gca,'FontSize',16); caxis([-0.1,0.1]);
title(sprintf('Preferred phase: T_{%g}-T_{%g}',Tarray(2),Tarray(1)));
hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
leftLim = [0.3,0.42,0.25,0.33];
axis(leftLim);

subplot(1,2,2);
patch(U_patch,V_patch,mydiffThresh,'LineStyle','None');
%patch(U_patch,V_patch,prefPhaseSmoCell{2}-prefPhaseSmoCell{1},'LineStyle','None');
axis equal;axis tight;colormap(hsvTrunc);
colorbar; set(gca,'FontSize',16); caxis([-0.1,0.1]);
title(sprintf('Preferred phase: T_{%g}-T_{%g}',Tarray(2),Tarray(1)));
%title('Preferred phase: T_{1.5}-T_3');
hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
rightLim = [0.34,0.47,0.65,0.75];
%rightLim = [0.33,0.48,0.74,0.84];
axis(rightLim)

pause;


%% NOT thresholded 
figure;subplot(1,2,1);
%patch(U_patch,V_patch,prefPhaseSmoCell{2}-prefPhaseSmoCell{1},'LineStyle','None');
patch(U_patch,V_patch,mydiff,'LineStyle','None');
axis equal;axis tight;colormap(hsvTrunc);
colorbar; set(gca,'FontSize',16); caxis([-0.1,0.1]);
title(sprintf('Preferred phase: T_{%g}-T_{%g}',Tarray(2),Tarray(1)));

hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
leftLim = [0.3,0.42,0.25,0.33];
axis(leftLim);

subplot(1,2,2);
patch(U_patch,V_patch,mydiff,'LineStyle','None');
%patch(U_patch,V_patch,prefPhaseSmoCell{2}-prefPhaseSmoCell{1},'LineStyle','None');
axis equal;axis tight;colormap(hsvTrunc);
colorbar; set(gca,'FontSize',16); caxis([-0.1,0.1]);
title(sprintf('Preferred phase: T_{%g}-T_{%g}',Tarray(2),Tarray(1)));

hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
rightLim = [0.34,0.47,0.65,0.75];
%rightLim = [0.33,0.48,0.74,0.84];
axis(rightLim)

pause;



%% Just check: in ortho vs. retro?
prefPhaseCycle = cell(size(prefPhaseSmoCell));
for j1=1:length(prefPhaseSmoCell)
    prefPhaseCycle{j1}= prefPhaseSmoCell{j1}>0.5;
end

%% Plot a binary switch
prefCycl=figure;
patch(U_patch,V_patch,prefPhaseCycle{2}-prefPhaseCycle{1},'LineStyle','None');
axis equal;axis tight;colormap(gray);
colorbar; set(gca,'FontSize',16); %caxis([-0.5,0.5]);
title(sprintf('Preferred phase: T_{%g}-T_{%g}',Tarray(1),Tarray(2)));
hold on;
line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2);

pause;

%% Plot point by point
prefComp=figure;
% Jitter randomly so we can see multiple points
jitt1 = prefPhaseSmoCell{1}+0.0005*randn(size(prefPhaseSmoCell{1}));
jitt2 = prefPhaseSmoCell{2}+0.0005*randn(size(prefPhaseSmoCell{1}));
plot(jitt1,jitt2,'.');
hold on;
plot(jitt1(Bdy2D.faceList),jitt2(Bdy2D.faceList),'.','MarkerSize',10);
plot([0 1],[0 1],'k--');
plot([0 1],[0.5 0.5],'m--');
plot([0.5 0.5],[0 1],'m--');
plot([0 1],[0.25 0.25],'c--');
plot([0.25 0.25],[0 1],'c--');
plot([0 1],[0.75 0.75],'r--');
plot([0.75 0.75],[0 1],'r--');
set(gca,'FontSize',20); axis equal;axis tight;
xlabel(sprintf('Preferred phase: T=%g s',Tarray(1)));
ylabel(sprintf('Preferred phase: T=%g s',Tarray(2)));
legend('All faces','Olfactory faces only');

Ortho_to_Retr = find(prefPhaseCycle{2}-prefPhaseCycle{1}==1);
Retr_to_Ortho = find(prefPhaseCycle{2}-prefPhaseCycle{1}==-1);


pause;

%% Identify certain sets of interest. We will then intersect w/ 
%% olfactory region

% Ortho peak for T=6, but retro peak for T=1.5
set1 = find(jitt1<0.3 & jitt2>0.5);

% Both ortho peak, but there is that slight delay
set2 = find(jitt1>0.2 & jitt1 <0.26 & jitt2>0.2 & jitt2 <0.26);


% Both are retro, and this is the "inertial delay" 
% I think I see
set3 = find(jitt1>0.72 & jitt2 > 0.72);

% 
set4 = find(jitt1>0.7 & jitt1 < 0.8 & jitt2 > 0.5 & jitt2 <= 0.72);

% Goes from Retro -> Ortho
set5 = find(jitt1>0.7 & jitt1 < 0.8 & jitt2 < 0.3);

maskSets = cell(5,1);
for j1=1:5
    maskSets{j1}=zeros(size(prefPhaseSmoCell{1}));
end
maskSets{1}(set1) = 1;
maskSets{2}(set2) = 1;
maskSets{3}(set3) = 1;
maskSets{4}(set4) = 1;
maskSets{5}(set5) = 1;

titleMaskSets=cell(5,1);
titleMaskSets{1}='Ortho peak for T=6, retro peak for T=1.5';
titleMaskSets{2}='Inertial delay during ortho phase';
titleMaskSets{3}='Inertial delay during retro phase';
titleMaskSets{4} = 'Slow at peak, but fast spread through early retro phase';
titleMaskSets{5}='Retro peak for T=6, ortho peak for T=1.5';

% Intersections
iset1 = intersect(set1,Bdy2D.faceList);
iset2 = intersect(set2,Bdy2D.faceList);
iset3 = intersect(set3,Bdy2D.faceList);
iset4 = intersect(set4,Bdy2D.faceList);
iset5 = intersect(set5,Bdy2D.faceList);

for j1=1:5
    figure; subplot(1,2,1);
    patch(U_patch,V_patch,1-maskSets{j1},'LineStyle','None');
    colormap(gray); hold on; axis equal;axis tight;
    set(gca,'FontSize',16);
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2); 
    axis(leftLim);
    title(titleMaskSets{j1});

    subplot(1,2,2);
    patch(U_patch,V_patch,1-maskSets{j1},'LineStyle','None');
    colormap(gray); hold on; axis equal;axis tight;
    set(gca,'FontSize',16);
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','g','LineWidth',2); 
    axis(rightLim);
end



if(0)

ifleft = 0;
ifright = 0;
if (ifleft)
    leftLim = [0.3,0.42,0.25,0.33];
    axis(leftLim);
elseif (ifright)
    rightLim = [0.34,0.47,0.65,0.75];
    %rightLim = [0.33,0.48,0.74,0.84];
    axis(rightLim)
end

end

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

if (0)
    %% 3D
figure;
patch(X_patch,Y_patch,Z_patch,prefPhaseSmoCell{1},'LineStyle','None');
axis equal;axis tight; colormap(RBmap);
colorbar; set(gca,'FontSize',16); clim([0,1]);
hold on; 
line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','g','LineWidth',2);
title('Pref phase: 3 s');
end

for j1=1:length(dirlist)
    %% plot this AGAIN
    figure;subplot(1,2,1);
    patch(U_patch,V_patch,prefPhaseSmoCell{j1},'LineStyle','None');
    axis equal;axis tight;colormap(RBmap);
    colorbar; set(gca,'FontSize',16); caxis([0,1]);
    title(sprintf('Preferred phase: T=%g',Tarray(j1)));
    hold on;
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
    leftLim = [0.3,0.42,0.25,0.33];
    axis(leftLim);
    
    subplot(1,2,2);
    patch(U_patch,V_patch,prefPhaseSmoCell{j1},'LineStyle','None');
    axis equal;axis tight;colormap(RBmap);
    colorbar; set(gca,'FontSize',16); caxis([0,1]);
    title(sprintf('Preferred phase: T=%g',Tarray(j1)));
    hold on;
    line(Bdy2D.X,Bdy2D.Y,Bdy2D.Z,'Color','k','LineWidth',2);
    rightLim = [0.34,0.47,0.65,0.75];
    %rightLim = [0.33,0.48,0.74,0.84];
    axis(rightLim)
end


figure;
patch(X_patch,Y_patch,Z_patch,prefPhaseSmoCell{2},'LineStyle','None');
axis equal;axis tight; colormap(RBmap);
colorbar; set(gca,'FontSize',16); clim([0,1]);
hold on; 
line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','g','LineWidth',2);
title('Pref phase: 1.5 s');


figure;
patch(X_patch,Y_patch,Z_patch,prefPhaseSmoCell{2}-prefPhaseSmoCell{1},'LineStyle','None');
axis equal;axis tight; colormap(RBmap);
colorbar; set(gca,'FontSize',16); clim([-0.1,0.1]);
hold on; 
line(Bdy3D.X,Bdy3D.Y,Bdy3D.Z,'Color','k','LineWidth',2);
title('Pref Phase: Diff');


if (0)
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

