function [w,wallInfo] = read_wss_mag_OF_face(filename,wallInfo)
% read_wss_mag_OF_face:  read an Open Foam wall shear stress (WSS) file
%       and return the MAGNITUDE of the WSS vector
%
% INPUT
%  filename :  wallShearStress file from Open Foam
% 
% OUTPUT
%  w  : nx3 matrix of WSS for each wall face
%  wallInfo  : supplementary info from file,
%      characterizing each segment of wall
%      (can be used later to validate match with location information)
%     wallInfo.Name:    name of boundary region
%     wallInfo.isWall:  1/0 flag to indicate wall/not wall
%     wallInfo.nFace:   number of faces in each region
%

fid = fopen(filename);
A = regexp(fileread(filename),'\n','split');

boundaryStart = find(contains(A,'boundaryField'),1);
bracelinenum  = find(contains(A,'{'));
bracelinenum  = bracelinenum(bracelinenum > boundaryStart+1);

% Examples of boundary regions we wish to save vs. not save
%
% SAVE
% boundaryField
% {
%     septal
%     {
%         type            calculated;
%         value           nonuniform List<vector> 
% 7808
% (
% (-8.28566e-06 2.99975e-06 -2.83666e-06)
%
%
% DO NOT SAVE
%     nostril
%     {
%         type            calculated;
%         value           uniform (0 0 0);
%     }


% Each is a boundary region
% Iterate through to check whether we need to read
%   values
nReg      = length(bracelinenum);
nameArray = cell(nReg,1);
isWall    = zeros(nReg,1);
nFace     = isWall;
for j=1:nReg
    lnum         = bracelinenum(j);
    nameArray{j} = strtrim(A{lnum-1});
    
    if find(regexp(A{lnum+2},'nonuniform'))
        isWall(j)=1;
        nFace(j) = str2num(A{lnum+3});
    end
    
end
  
% Now, iterate through wall regions and read WSS vectors
allWSS = zeros(sum(nFace),3);
startI = 1;
for j=1:nReg
    if (isWall(j))   %Otherwise skip
        % Location in array to place points
        endI = startI + nFace(j)-1;
        
        % Find line in file where we will start to read
        % See example above
        startLine = bracelinenum(j)+5;
        for j1=1:nFace(j)
            allVal = textscan(A{j1+startLine-1},'(%f %f %f)');
            allWSS(startI+j1-1,:) = [allVal{1} allVal{2} allVal{3}];
        end
        
        % Reset startI for the next region
        startI = startI + nFace(j);
    end
end

% % CHECK
% startI = 1;
% for j=1:nReg
%     if (isWall(j))
%         endI = startI+nFace(j)-1;
% %        [allWSS(startI,:);allWSS(endI,:)]
% %        [startI endI]
%         
%         startI = startI+nFace(j);
%     end
%         
% end

wallInfo.Name       = nameArray;
wallInfo.isWall     = isWall;
wallInfo.nFace      = nFace;
% disp(nameArray{1})
% isWall
% nFace

% Return magnitude
w=  sqrt(sum(allWSS.^2,2));


end