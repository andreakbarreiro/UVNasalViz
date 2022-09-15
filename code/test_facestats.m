
% test_facestats.m
addpath('util');

caseDir = '../Elad_Test_Case/';

% Boundary file name
fname   = 'constant/polymesh/boundary';

[fp fstats]=facestats([caseDir fname]);
