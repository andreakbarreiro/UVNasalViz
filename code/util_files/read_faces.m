function faces = read_faces(filename)
%% read_faces: read faces from an Open Foam file
%
% % !!!Note:  All faces MUST be triangles
%
% INPUTS
%  filename :       must be a "faces" file produced by Open Foam
%                   (usually found in constant/polymesh/)
% OUTPUTS
%  faces :        (n x 3) array of corners for each face
%

% We read in each line
A = regexp(fileread(filename),'\n','split');

% Find first line with an "("
linenum = find(contains(A,'('),1);

% The previous line contains the number of faces
nFaces        = str2num(A{linenum-1});

% line on which faces will start
startLine   = linenum+1;

faces = zeros(nFaces, 3);
for j=1:nFaces
    allCoords=textscan(A{j+startLine-1},'3(%d %d %d)');
    faces(j,:) = [allCoords{1} allCoords{2} allCoords{3}];
end

% % From Abdullah
% if (0)
% facefile = filename;
% 
% opts = detectImportOptions(facefile);
% 
% opts.SelectedVariableNames = {'Var1', 'Var2', 'Var3'};
% 
% T = readtable(facefile,opts);
% 
% T(end-1:end,:) = [];
% 
% T.Var1 = erase(T.Var1,'3(');
% T.Var3 = erase(T.Var3,')');
% 
% T.Var1 = str2double(T.Var1);
% T.Var3 = str2double(T.Var3);
% 
% 
% f = table2array(T);
% end

end