% get_sim_metadata.m
function allInfo = get_sim_metadata(SimName,metadata_dir)

% Must pass in metadata_dir.
% 
%metadata_dir = '/users/abarreiro/nasoCFD/';

metadir_fname = 'Index_InthaDomainSims.xlsx';

% Read file
% ',': CSV
% 1: number of header lines
A= readtable([metadata_dir metadir_fname],"FileType","spreadsheet",...
    "ReadRowNames",true, "ReadVariableNames",true,...
    "PreserveVariableNames", true,...
    "ImportErrorRule","fill");


% Which column is which? 
simName_colName = 'SimName';
path_colName = 'Path on M3';
freq_colName = 'Forcing Frequency (1/s)';
mesh_colName = 'Mesh';
time_colName = 'Total time (s)';
dt_colName   = 'dt (s)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Which simulation?

%SimName = 'intha_0.5mil_tetra';

% Instead of hard-coding file locations, period, etc.
% Use this lookup table
mainCaseDir = '/projects/abarreiro/nasoCFD/nasoCFD/'

% Where results are stored, relative to main case directory 
caseDir  = [mainCaseDir A{SimName,path_colName}{1} SimName];

if (caseDir(end) ~= '/')
    caseDir = [caseDir '/'];
end

% Where OpenFoam mesh is
% Boundary file name
bfname   = 'constant/polyMesh/boundary';

% Full path is: mainCaseDir + caseDir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Where to find UV file, and boundary file?

whichMesh=A{SimName,mesh_colName}{1};

%%  Assume the mesh ID determines the location
%   UV file will be in the same place as the boundary file
%   
topBdydir = '/users/abarreiro/nasoCFD/UVNasalViz-main/';
if (strcmp(whichMesh,'old_mesh'))
    bdydir  = [topBdydir 'caseDirs/Inthavong_Test_Case/'];
    objdir = bdydir; 
    objfile  = 'uv_FIXED.obj';
elseif (strcmp(whichMesh,'new_mesh'))
    % Attempted to improve mesh, in 2024
    bdydir  = [topBdydir 'caseDirs/Intha_newMesh/'];
    objdir = bdydir;
    objfile = 'uv_intha_0.5mil_fixU.obj';
elseif (strcmp(whichMesh,'fluent_0.5mil'))
    % Used Ansys, in 2025. Includes boundary layers
    bdydir  = [topBdydir 'caseDirs/Intha_BoundaryLayers/'];
    objdir = bdydir;
    objfile = 'uv_intha_0.5mil_tetra.obj';
end

bdyfile = 'bdy_All_Inthavong.mat';


allInfo = struct;
allInfo.caseDir = caseDir;

allInfo.bdydir  = bdydir;
allInfo.bdyfile = bdyfile;
allInfo.objdir = objdir;
allInfo.objfile = objfile;

%% Info that is specific to the simulation, not the mesh
%   These are numbers, they will not be contained in a cell
allInfo.freq  = A{SimName,freq_colName};
allInfo.lastT = A{SimName,time_colName};
allInfo.dt    = A{SimName,dt_colName};