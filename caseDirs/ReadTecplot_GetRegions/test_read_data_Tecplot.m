
% Tecplot.360 file
%
% Subdivided by regions
tecfile = 'DUN001_mapping_example_tec.dat';

% Comparison .obj file that has exterior points, but not labeled with
% regions
if (0)
objdir  = '../../UVViz/caseDirs/Inthavong/';
%objfile = 'Geometry_mapped_d4_Closed.obj';
objfile = 'Geometry_mapped_d4.obj';
end

objdir  = '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/';
objfile = 'uv2.obj';

[ROld,P,FACE_Old,S,lmax,zo_Info]=read_data_Tecplot_wZones(tecfile);
figure;
subplot(1,2,1);
plot3(ROld(:,1),ROld(:,2),ROld(:,3),'.');
set(gca,'FontSize',14);
title('From tecplot file')

% Read in obj file
[r,face,r2R,R]=read_3D_from_OBJ([objdir objfile]);

subplot(1,2,2);
plot3(R(:,1),R(:,2),R(:,3),'.');
set(gca,'FontSize',14)
title('From .obj file')


if (1)
% Here is a curious thing. There appear to be all kinds of unused points.
% to visualize, look at the following:
% 
maxPts = max(r2R);

% Finally, the ACTUALLY USED points
usedpts = unique(r2R);
unusedpts = setdiff(1:maxPts,usedpts); 

figure;
plot3(R(usedpts,1),R(usedpts,2),R(usedpts,3),'.');

%plot3(R(1:maxPts,1),R(1:maxPts,2),R(1:maxPts,3),'.');

hold on;
pause;

% They appear to be a TRANSLATION of the actual points.
plot3(R(maxPts+1:end,1),R(maxPts+1:end,2),R(maxPts+1:end,3),'.');

pause;

plot3(R(unusedpts,1),R(unusedpts,2),R(unusedpts,3),'.');
set(gca,'FontSize',14);
title('From .obj file, in 3 sets')
legend('Used','Extra','Unused')
end


fprintf('length of r = %d\n ', length(r))
fprintf('length of R = %d\n', length(R))
fprintf('length of r2R = %d\n ', length(r2R))
fprintf('length of face = %d\n', length(face))

% 
if (1)
% Remove any points that don't map to anything
maxPts = max(r2R); 
R = R(1:maxPts,:);
end


% NOW plot them together

figure;
plot3(ROld(:,1),ROld(:,2),ROld(:,3),'.');
set(gca,'FontSize',14);
hold on;
plot3(R(usedpts,1),R(usedpts,2),R(usedpts,3),'.');
plot3(R(unusedpts,1),R(unusedpts,2),R(unusedpts,3),'.');
legend('Tecplot','OBJ','OBJ, not used');

disp('Min/Max for ROld')
[min(ROld);max(ROld)]
diff([min(ROld);max(ROld)])

mROld=mean([min(ROld);max(ROld)])

disp('Min/Max for R')
[min(R(usedpts,:));max(R(usedpts,:))]
diff([min(R(usedpts,:));max(R(usedpts,:))])

mR=mean([min(R(usedpts,:));max(R(usedpts,:))])

toAdd = mR-mROld;

%% Now translate mROld by the needed amount
ROld = ROld + toAdd;


disp('NEW Min/Max for ROld')
[min(ROld);max(ROld)]

disp('NEW Min/Max for R')
[min(R(usedpts,:));max(R(usedpts,:))]

%% Now, assign a zone to "faces"
nFaces = length(face);

centroids = zeros(nFaces,3);
for k=1:nFaces
    % Find the centroid of the face in 3D
    centroids(k,:)=mean(R(r2R(face(k,:)),:));
end
zMatchFac = zeros(nFaces,3);

%% Since TECplot file is organized by zones, use ROld index
%% to assign zone
firstInds = cumsum(zo_Info{2})+1;
firstInds = [1;firstInds];

for k=1:nFaces
    currP = centroids(k,:); %currP

    % Find the closest point(s) in ROld
    distP = sum((ROld-currP).^2,2);

    myind = find(distP == min(distP));
    
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

% for k=1:nFaces
%     zMatchFac(k,:) = sort(reshape(zoneMatch(face(k,:),:),1,9),'descend');
% end

%% If any of the closest points is in either of the olfactory 
%% regions
touchFacInd = (any(zMatchFac==7,2)) | (any(zMatchFac==8,2));
whichFac = find(touchFacInd);

%% Find boundary: both in 3D, and 2D

[lb3,X3,Y3,Z3] = findb(R,r2R(face(whichFac,:)));

[lb2,X2,Y2,Z2]= findb(r,face(whichFac,:));



addpath('../../UVViz/code/util_files');
addpath('../../UVViz/code/util_geom');

% Where to find a WSS file
caseDir = '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/'
timestep    = 0.5;
wssfname    = sprintf('%s%g/wallShearStress',caseDir,timestep);

w = read_wss_mag_OF_face(wssfname);

% w should be defined ON THE SAME FACES
% i.e. there is a 1-1 mapping between faces in the 3D model
% and in the UV-mapping.

%% 3D figure
figure;
x=R(:,1); y=R(:,2);  z=R(:,3);

FACE = r2R(face);
X_patch=x(FACE)';Y_patch=y(FACE)'; Z_patch=z(FACE)';
patch(X_patch,Y_patch,Z_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;
 
hold on;
line(X3,Y3,Z3,'Color','r','LineWidth',2);

%% 2D figure
figure;
u=r(:,1); v=r(:,2); 
U_patch=u(face)';V_patch=v(face)';
patch(U_patch,V_patch,w','LineStyle','None');
axis equal;axis tight;
colorbar;

hold on;
line(X2,Y2,Z2,'Color','r','LineWidth',2);


%% Test boundary
% Make an indicator function for region
wIndic = zeros(size(w));
wIndic(whichFac) = 1;

figure;
patch(X_patch,Y_patch,Z_patch,wIndic','LineStyle','None');
axis equal;axis tight;
colorbar;
 
hold on;
line(X3,Y3,Z3,'Color','r','LineWidth',2);

figure;
patch(U_patch,V_patch,wIndic','LineStyle','None');
axis equal;axis tight;
colorbar;

hold on;
line(X2,Y2,Z2,'Color','r','LineWidth',2);


if (1)
    Bdy3D = struct(); Bdy3D.lb = lb3;
    Bdy3D.X = X3; Bdy3D.Y = Y3; Bdy3D.Z = Z3;

    % This should include all faces inside
    Bdy3D.faceList =  whichFac;

    Bdy2D = struct(); Bdy2D.lb = lb2;
    Bdy2D.X = X2; Bdy2D.Y = Y2; Bdy2D.Z = Z2; 
    
    Bdy2D.faceList =  whichFac;

    % Now save
    sourcefile = [objdir objfile];
    bdyfile    = 'bdy_Olfac_Inthavong.mat';
    save(bdyfile,'sourcefile','Bdy3D','Bdy2D','toAdd');
end