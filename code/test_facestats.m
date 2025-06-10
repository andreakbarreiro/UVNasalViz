
% test_facestats.m
addpath('util_files');

caseDir = '../caseDirs/Elad_Test_Case/';

% Boundary file name
fname   = 'constant/polyMesh/boundary';

[fp fstats]=facestats([caseDir fname]);
