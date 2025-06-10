% Testing on M3
%   Feb 2025
%


addpath('util_files');
addpath('util_geom');
addpath('util_plot');

%% Which mesh? Old or new
%% Old= circa 2022, used to do first PPMs in 2023
%% New=circa 2024, supposed to be corrected for nonorthogonality
use_new_mesh_flag = 0;

% Read both geometries

% ON M3
caseDir_Head = '/projects/abarreiro/nasoCFD/nasoCFD/intha/'
caseDir_List = {[caseDir_Head 'new_mesh/intha_0.5mil_fixU/'];...
        [caseDir_Head 'old_mesh/intha_normal/']};

bdydir_List = {'../caseDirs/Inthavong_newMesh/';...
        '../caseDirs/Inthavong_Test_Case/'};
bdyfile_List = {'bdy_All_Inthavong_TEST.mat';...
        'bdy_All_Inthavong.mat'};

% These come from OF mesh
R_OF    = cell(length(bdyfile_List),1);
FACE    = cell(length(bdyfile_List),1);
fp      = R_OF;
fstats  = R_OF;

% These come from UV file
R       = cell(length(bdyfile_List),1);
r       = R;
face    = R;
r2R     = R;

Bdy3D_All = R;
Bdy2D_All = R;

for j1=1:length(bdydir_List)
    % Read OF points
    fname   = 'constant/polyMesh/boundary';
    [fp{j1} fstats{j1}]=facestats([caseDir_List{j1} fname]);
    [R_OF{j1} FACE{j1} S lmax] = read_polyMesh(caseDir_List{j1});

    % Read UV file
    load([bdydir_List{j1} bdyfile_List{j1}]);

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

    Bdy3D_All{j1}=Bdy3D; Bdy2D_All{j1}=Bdy2D;
    %% "sourcefile" will contain the proper UV file name
    %%    However, directory is probably different from where file was produced
    %%    
    %% Strip directory off
    uvFname = split(sourcefile,'/'); uvFname = uvFname{end};
    % Assume located in same place as boundary file
    uvDir = bdydir_List{j1};
   
    [r{j1} face{j1} r2R{j1} R{j1}]=read_3D_from_OBJ([uvDir uvFname]);


end

 
for j1=1:length(caseDir_List)
    figure;
    plot3(R_OF{j1}(:,1),R_OF{j1}(:,2), R_OF{j1}(:,3),'.'); hold on;
    plot3(R{j1}(:,1),R{j1}(:,2), R{j1}(:,3),'.'); 
    set(gca,'FontSize',16);
    legend('Open Foam','UV');
    title(sprintf('Mesh \# %d, j1'));
end

figure;
for j1=1:2
    plot3(R_OF{j1}(:,1),R_OF{j1}(:,2), R_OF{j1}(:,3),'.'); hold on;
end
set(gca,'FontSize',16);
legend('Open Foam #1','Open Foam #2');

figure;
for j1=1:2
    plot3(R{j1}(:,1),R{j1}(:,2), R{j1}(:,3),'.'); hold on;
end
set(gca,'FontSize',16);
legend('UV \#1','UV \#2');

for j1=1:2
    if (j1==1)
        disp('New mesh')
    else
        disp('Old mesh')
    end
    disp('Size of R vs R_OF');
    disp([length(R{j1}) length(R_OF{j1})])

    disp('Size of r vs r2R');
    disp([length(r{j1}) length(r2R{j1})])

    disp('Unique points in r2R')
    disp(length(unique(r2R{j1})))

    disp('faces vs. FACES')
    disp([length(face{j1}) length(FACE{j1})])

    disp('facestats')
    disp(fp{j1})
end

% (1) Create modified "R" arrays that do not have open boundaries
Rtemp = cell(size(R));
targetCell = Rtemp;
for j1=1:2
    target_array = sort(unique(r2R{j1}));
    Rtemp{j1} = R{j1}(target_array,:);
    % These are the points that were maintained.
    % Must reference back to these
    targetCell{j1} = target_array;
end
ptCloudA = pointCloud(Rtemp{1}); % fixed
ptCloudB = pointCloud(Rtemp{2}); % moving


[tform,movingReg] = pcregistericp(ptCloudB,ptCloudA,Metric="planeToPlane");
%ROld = movingReg.Location;

newptCloudB = movingReg.Location;

%(2) For each pointB, find nearest point in cloud A
nptsA = length(ptCloudA.Location);
nptsB = length(newptCloudB);
closestAtoB = zeros(nptsB,1);
for j1=1:nptsB
    allDist2 = sum((ptCloudA.Location-repmat(newptCloudB(j1,:),nptsA,1)).^2,2);
    closestAtoB(j1) = find(allDist2 == min(allDist2));
end

temp = [newptCloudB ptCloudA.Location(closestAtoB,:)];

% (3) Now find matching UV point(s)
r2Rtemp=cell(size(r2R));
for j1=1:2
    % Add on "r" index
    r2Rtemp{j1} = [[1:length(r2R{j1})]' r2R{j1}];
end
R2rtemp = cell(size(r2R));
for j1=1:2
    R2rtemp{j1} = sortrows(r2Rtemp{j1},2);
    R2rtemp{j1}(:,2)
end

uvA = zeros(length(newptCloudB),2);
uvB = uvA;

% Make sure this reflects the right R2r
R2rB = R2rtemp{2}; rB = r{2};
R2rA = R2rtemp{1}; rA = r{1};

targetB  = targetCell{2}; targetA = targetCell{1};

% (4) Find UV coordinates. For each "pointB",
% find (4.1) corresponding UV point in 2D
%  and (4.2) corresponding UV point of closest "pointA"
for j1=1:length(newptCloudB)
    % Look for where the 'R' index is in column 2
    %   Column 1 should give the 'r' index
    % 
    % To index into r2R; make sure to use 
    % targetCell
    whichInd  = find(R2rB(:,2)==targetB(j1),1);
    if (isempty(whichInd))
        uvB(j1,:) =[nan nan];
    else
        uvB(j1,:) = rB(R2rB(whichInd,1),:);
    end
    % First, find index of closest "A" point
    whichInd  = find(R2rA(:,2)==targetA(closestAtoB(j1)),1);
    if (isempty(whichInd))
        uvA(j1,:) =[nan nan];
    else
        uvA(j1,:) = rA(R2rA(whichInd,1),:);
    end
end

% Finally remove nan rows
anyNan = find(any(isnan([uvA uvB]),2));

uvA(anyNan,:)=[];
uvB(anyNan,:)=[];

figure;
subplot(2,1,1);
plot(uvA(:,1),uvB(:,1),'.');
subplot(2,1,2);plot(uvA(:,2),uvB(:,2),'.');


% Look at a square
uLim = [0.3 0.35]; vLim=[0.37 0.4];
temp = find(uvA(:,1)>=uLim(1)&uvA(:,1)<=uLim(2) ...
    & uvA(:,2)>=vLim(1)&uvA(:,2)<=vLim(2));
figure;
subplot(1,2,1);
plot(uvA(temp,1),uvA(temp,2),'.'); hold on;
plot(uvB(temp,1),uvB(temp,2),'.');


uLimB = [min(uvA(temp,1)) max(uvA(temp,1))];
vLimB = [min(uvA(temp,2)) max(uvA(temp,2))];

uptsSquare = [uLimB(1), vLimB(1);...
    uLimB(1), vLimB(2);...
    uLimB(2),vLimB(2);...
    uLimB(2), vLimB(1);...
    uLimB(1), vLimB(1)];

plot(uptsSquare(:,1),uptsSquare(:,2));

uvNewMesh = uvA;
uvOldMesh = uvB;

save('uvMatch_Old_vs_New.mat','uvNewMesh','uvOldMesh');