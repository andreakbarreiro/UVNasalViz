% Compare OBJ files.
%
% I have TWO of the same region
%  Unsure of their provenance
%   UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/uv2.obj
%   Timestamped: Oct 22 2022
%  (This was not long after I was creating the mesh for this region:
%     Oct 16 2022)
%  
%   These files refer to "Geometry_mapped_d4_Closed.obj"
%     which is now in
%       UVViz/caseDirs/Inthavong
%   It appears to be a scaled down (smaller) version
%      of Geometry_mapped_Orig.obj (dated July 2021, 10.2 MB)
%
%  UVViz/caseDirs/intha_amp_50/uv.obj
%  Timstamped: Jan 26 2023
%
% 
% UPDATE: After investigation
%   uv.obj is a perfect translation of the OpenFoam points
%      (modulo the "extra copy" of the boundary points)
%
%   uv2.obj has a different translation for wall vs. patch 
%      regions (nostrils, nasopharynx)
% The numerical values of the translations do not appear to be related to
% one another.
%
%
% QUESTION: HOW DO I FIX THIS WITH THE LEAST AMOUNT OF WORK
%   POSSIBLE?????
%  Here is what I will do.
%  - Write out new OBJ file.
%  -  AFAIK both uv and uv2 contain valid unwrappings, 
%       albeit with REVERSED ORIENTATIONS (i.e. left-right, I think, have
%       been swapped).
%  - FOR BOTH: remove extra points
%  - FOR BOTH: find translation relative to OpenFoam points 
%  -
%   
%
addpath('../../UVViz/code/util_files');
addpath('../../UVViz/code/util_geom');
addpath('../../UVViz/code/util_plot');


% Let us try to revisit these few files
fdirlist = {'../../UVViz/caseDirs/Inthavong/';...
    '../../UVViz/caseDirs/Inthavong/';...
    '../../UVViz/caseDirs/intha_amp_50/';...
    '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/';...    
    };
fnamelist = {'Geometry_mapped_Orig.obj';...
    'Geometry_mapped_d4_Closed.obj';...
    'uv.obj';...
    'uv2.obj';...
    };

nFiles = length(fdirlist);

r = cell(nFiles,1);
R = r;
r2R = r;
face = r;

for j=1:nFiles
    % Read in obj file
    [r{j},face{j},r2R{j},R{j}]=read_3D_from_OBJ([fdirlist{j} fnamelist{j}]);
end

% Get rid of "extra points"
for j=1:nFiles
    maxRind = max(r2R{j});
    R{j}(maxRind+1:end,:)=[];
end


%% Read OF geometry, if not already done
read_polyMesh_flag = 1;
%read_uv_flag       = 0;
if (read_polyMesh_flag)
    % I don't love having to keep multiple copies of the polyMesh
    % and UV files. However, it's probaly better than having endless
    % confusion
    %
    % Inthavong directory
    caseDir = '../../UVViz/caseDirs/intha_amp_50/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp,fstats]=facestats([caseDir fname]);
    [R_of,FACE_of,S_of,lmax_of] = read_polyMesh(caseDir);
end


figure; 
for j=1:nFiles
    subplot(2,2,j);
    plot3(R{j}(:,1),R{j}(:,2),R{j}(:,3),'.');
    set(gca,'FontSize',14);
end

figure;
for j=1:nFiles
   subplot(2,2,j);
   plot(r{j}(:,1),r{j}(:,2),'.');
   set(gca,'FontSize',14);
end



 temp=R{3}-R{4};

 maxPts = length(temp);

 temp = round(temp,4);
 set1 = find(temp(:,1)==temp(1,1));
 set2 = setdiff([1:maxPts]',set1);  

 tempB = R{3}-R_of;

 maxPtsB = length(tempB);

 tempB = round(tempB,4);
 set1B = find(tempB(:,1)==tempB(1,1));
 set2B = setdiff([1:maxPtsB]',set1B);


 figure;
 subplot(1,3,1); plot3(R{3}(set1,1),R{3}(set1,2),R{3}(set1,3),'.');
 hold on;plot3(R{3}(set2,1),R{3}(set2,2),R{3}(set2,3),'.');

 subplot(1,3,2);
 plot3(R{4}(set1,1),R{4}(set1,2),R{4}(set1,3),'.');
 hold on;plot3(R{4}(set2,1),R{4}(set2,2),R{4}(set2,3),'.');