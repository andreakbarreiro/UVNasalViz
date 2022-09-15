function pcoords = read_points(filename)
%% read_points: read point coordinates from an Open Foam file 
%
% INPUTS
%  filename :       must be a "points" file produced by Open Foam
%                   (usually found in constant/polymesh/)
% OUTPUTS
%  pcoords :        (n x 3) array of 3D coordinates

% We read in each line
A = regexp(fileread(filename),'\n','split');

% Find first line with an "("
linenum = find(contains(A,'('),1);

% The previous line contains the number of points
npts        = str2num(A{linenum-1});

% line on which points will start
startLine   = linenum+1;

pcoords = zeros(npts, 3);
for j=1:npts
    allCoords=textscan(A{j+startLine-1},'(%f %f %f)');
    pcoords(j,:) = [allCoords{1} allCoords{2} allCoords{3}];
end



% % Code from Abdullah
% if (0)
% pointfile = filename;
% 
% opts = detectImportOptions(pointfile);
% 
% opts.SelectedVariableNames = {'Var1', 'Var2', 'Var3'};
% 
% T = readtable(pointfile,opts);
% 
% T(end-1:end,:) = [];
% 
% T.Var1 = erase(T.Var1,'(');
% T.Var3 = erase(T.Var3,')');
% 
% T.Var1 = str2double(T.Var1);
% T.Var3 = str2double(T.Var3);
% 
% 
% p = table2array(T);
% 
% end
end