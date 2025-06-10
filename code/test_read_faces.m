%test_read_faces
%
% Test "read_faces"
%In particular: new changes to read 
% %  in 4-sided (or greater) faces

addpath('util_files');
addpath('util_geom');
addpath('util_plot');

if (0)
% Old mesh: only 3 sided faces
caseDir =  '/projects/abarreiro/nasoCFD/nasoCFD/intha/old_mesh/intha_normal/'

else
% New mesh w/ boundary layers
% Some faces have 4 sides
caseDir =  '/projects/abarreiro/nasoCFD/nasoCFD/test_sims/intha_0.5mil_tetra/'
end

[fp fstats] = facestats([caseDir 'constant/polyMesh/boundary']);

% For each physical region
startFaceIndex  =  fp(:,1) + 1;     %Offset by 1 because OpenFoam 
                                        %starts lists w/ "0"
nFace           =  fp(:,2);

% Read faces
% 
faces = read_faces([caseDir 'constant/polyMesh/faces']);
    
% Determine which boundary pieces are WALLS
whichWalls = find(fp(:,3));
     