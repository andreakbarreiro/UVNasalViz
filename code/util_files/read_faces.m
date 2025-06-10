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

% Make larger than needed. Remove unneeded columns at the end
maxNumCor = 6;

% Use nan
faces = nan(nFaces, maxNumCor);

for j=1:nFaces
    % Each line begines with an integer, This is the # of corners for the
    % face
    %disp(A{j+startLine-1})

    % For some reason, this returns 2 integers:
    %   the number of corners AND the first index
    nPts  = textscan(A{j+startLine-1},'%d( ');
    % To get only the number of corners
    nPts = nPts{1}; nPts=nPts(1);
    if (nPts > maxNumCor)
        print('Warning! You need to increase size of "faces" ');
    else
        % Will return an array of cells
        % Each cell contains an array with a single integer (which is the
        % point index)
        allCoords=textscan(A{j+startLine-1}, ...
            [num2str(nPts) '(' repmat('%d ',1,nPts-1) '%d)']);
      
    end
  
    % Will turn into a single row of nPts integers
    faces(j,1:nPts) = cell2mat(allCoords);
    %[allCoords{1} allCoords{2} allCoords{3}];
end

% Remove unnecessary columns.

toremove=find(all(isnan(faces)));
faces(:,toremove)=[];

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