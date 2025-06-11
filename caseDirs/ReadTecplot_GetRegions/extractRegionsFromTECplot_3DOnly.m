% Associate faces to regions from
% TECPLOT file (provided by Inthavong, etc..
%

% So we can access "get_sim_metadata.m"
% Need to iscover the full path vs. hard-coding it here
addpath('/users/abarreiro/nasoCFD/UVNasalViz-main/code');
addpath('/users/abarreiro/nasoCFD/UVNasalViz-main/code/util_files');
addpath('/users/abarreiro/nasoCFD/UVNasalViz-main/code/util_geom');



% Subdivided by regions
tecfile = 'DUN001_mapping_example_tec.dat';

% _3DOnly
%
% Issue is that OpenFoam mesh is not 1-1 with 
%  .obj file, depending on how the .obj file was created.
%  So, try to adapt to the OpenFoam file, so we can have a
%  version that works for 3D plots.
%
% I don't know why, AS creates the .obj files, and there seems to be no
% rhyme or reason to them. So we are forced to adapt
% 


use_simReg_flag = 1;

if (use_simReg_flag)
    % If you want to get file locations, etc.. from the "simulation
    % registration" file, an Excel spreadsheet
    metadata_dir = '/users/abarreiro/nasoCFD/';

    %% First thing: Identify which simulation
    SimName = 'intha_0.5mil_tetra';

    %% Get info about file locations, etc..
    allInfo = get_sim_metadata(SimName,metadata_dir);

    %% Please note: if bdyfile exists, we will rewrite it!
    % Make sure to save if this is not the wanted behavior
    warning('Boundary file will be rewritten! Make sure to save the file if this is not the desired behavior.')
else    
    caseDir_Head = '/projects/abarreiro/nasoCFD/nasoCFD/intha/'
    caseDir = [caseDir_Head 'new_mesh/intha_0.5mil_fixU/'];

    objdir = caseDir;
    objfile = 'uv_intha_0.5mil_fixU.obj'
    bdyfile    = 'bdy_All_Inthavong_TEST.mat';
   %objdir  = '../../UVViz/caseDirs/intha_amp_50/';
   %objfile = 'uv_FIXED.obj'; 
   allInfo = struct('caseDir',caseDir,...
        'bdydir', bdydir,...
        'bdyfile', bdyfile);
end

[ROld,P,FACE_Old,S,lmax,zo_Info]=read_data_Tecplot_wZones(tecfile);


% Read in obj file
%[r,face,r2R,R]=read_3D_from_OBJ([allInfo.objdir allInfo.objfile]);
[R,FACE,~,~] = read_polyMesh(allInfo.caseDir);

%fprintf('length of r = %d\n ', length(r))
fprintf('length of R = %d\n', length(R))
%fprintf('length of r2R = %d\n ', length(r2R))
fprintf('length of face = %d\n', length(FACE))


disp('Min/Max for ROld')
[min(ROld);max(ROld)]

disp('Min/Max for R')
[min(R);max(R)]

%% I may need to rotate points. Perform a "registration" 
%% so points align optimally.
% Create point cloud objects from your arrays

if (1)
    % Treat ROld as the moving points. We want our mesh
    % (i.e. the mesh wrt which our computation was performed) to be 
    % preserved. In contrast, this is the ONLY time we will refer to the 
    % coordinates in the tecplot file.
    % 
    % Create point cloud objects from your arrays
    ptCloudA = pointCloud(R); % fixed
    ptCloudB = pointCloud(ROld); % moving


    [~,movingReg] = pcregistericp(ptCloudB,ptCloudA,Metric="planeToPlane");
    ROld = movingReg.Location;
else
    % In Abdullah's scripts
    ptCloudA = pointCloud(ROld); % fixed
    ptCloudB = pointCloud(R); % moving

    [~,movingReg] = pcregistericp(ptCloudB,ptCloudA,Metric="planeToPlane");
    R = movingReg.Location;
end

disp('Min/Max for ROld')
[min(ROld);max(ROld)]

disp('Min/Max for R')
[min(R);max(R)]

%% Now, assign a zone to "faces"
nFaces = length(FACE);

centroids = zeros(nFaces,3);
for k=1:nFaces
    % Find the centroid of the face in 3D
    %centroids(k,:)=mean(R(r2R(face(k,:)),:));
    centroids(k,:)=mean(R(FACE(k,:),:));
end

zMatchFac = zeros(nFaces,3);
actMatchFac = zMatchFac;

%% Since TECplot file is organized by zones, use ROld index
%% to assign zone
firstInds = cumsum(zo_Info{2})+1;
firstInds = [1;firstInds];

for k=1:nFaces
    currP = centroids(k,:); %currP

    % Find the closest point(s) in ROld
    distP = sum((ROld-currP).^2,2);

    myind = find(distP == min(distP));
    actMatchFac(k,1:length(myind))=myind';

    %% Determine which regions those are in
    %% Could be up to 3 points
    zMatchFac(k,1)= find(firstInds>myind(1),1)-1;
    if (length(myind)>=2)
        zMatchFac(k,2)= find(firstInds>myind(2),1)-1;
    end
    if (length(myind)>=3)
        zMatchFac(k,3)= find(firstInds>myind(3),1)-1;
    end
end

%% FINAL SET
% First, default is that we pick first on list.
% If there is actually more than one: use mode
zMatchFinal = zMatchFac(:,1);


% Hypothesis: if there is more than 1 match,
% that is a point that lies on the bdy of two regions
%
sum((zMatchFac(:,2)>0))

sum((zMatchFac(:,3)>0))

actMatchFac(:,4:end)=[];

which3 = find(zMatchFac(:,3)>0);
% for j1=1:length(which3)
%     actMatchFac(which3(j1),:)
%     ROld(actMatchFac(which3(j1),:),:)
% end

% For each of these faces, get ALL points
% assoc. with face
for j1=1:length(which3)
    pts = FACE(which3(j1),:);
    allClose = [];
    for k1=1:3
        %currP = R(r2R(pts(k1)),:);
        currP = R(pts(k1),:);

        % Find the closest point(s) in ROld
        distP = sum((ROld-currP).^2,2);
        myind = find(distP == min(distP));
        allClose = [allClose;myind];
    end
    allCloseReg = zeros(size(allClose));
    for k1=1:length(allCloseReg)
        allCloseReg(k1)=find(firstInds>allClose(k1),1)-1;
    end
    %allCloseReg

    % Let's TRY the mode
    zMatchFinal(which3(j1))=mode(allCloseReg);
    %zMatchFac(which3(j1),:)
end

% Do something similar for "twos"
which2 = find(zMatchFac(:,2)>0 & ~zMatchFac(:,3));

for j1=1:length(which2)
    pts = FACE(which2(j1),:);
    allClose = [];
    for k1=1:3
        %currP = R(r2R(pts(k1)),:);
        currP = R(pts(k1),:);

        % Find the closest point(s) in ROld
        distP = sum((ROld-currP).^2,2);
        myind = find(distP == min(distP));
        allClose = [allClose;myind];
    end
    allCloseReg = zeros(size(allClose));
    for k1=1:length(allCloseReg)
        allCloseReg(k1)=find(firstInds>allClose(k1),1)-1;
    end
    allCloseReg

    % Let's TRY the mode
    zMatchFinal(which2(j1))=mode(allCloseReg);
    zMatchFac(which2(j1),:)
end

% Now get boundary for each region
nReg = length(zo_Info{2});

whichFac = cell(nReg,1);
for j1=1:nReg
    whichFac{j1}=find(zMatchFinal==j1);
end

allBdy3D = cell(nReg,1);
%allBdy2D = cell(nReg,1);
for j1=1:nReg
    [lb3,X3,Y3,Z3] = findb(R,FACE(whichFac{j1},:));
    %[lb2,X2,Y2,Z2] = findb(r,face(whichFac{j1},:));
    Bdy3D = struct(); Bdy3D.lb = lb3;
    Bdy3D.X = X3; Bdy3D.Y = Y3; Bdy3D.Z = Z3;
    
    % This should include all faces inside
    Bdy3D.faceList =  whichFac{j1};

    %Bdy2D = struct(); Bdy2D.lb = lb2;
    %Bdy2D.X = X2; Bdy2D.Y = Y2; Bdy2D.Z = Z2; 
    
    %Bdy2D.faceList =  whichFac{j1};

    allBdy3D{j1}=Bdy3D;
    %allBdy2D{j1}=Bdy2D;
end

% Now save

sourcefile = [allInfo.objdir allInfo.objfile];
save([allInfo.bdydir allInfo.bdyfile(1:end-4) '_3DOnly.mat'],'sourcefile',...
    'allBdy3D','zo_Info');


