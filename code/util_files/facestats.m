function [fp, fstats] = facestats(filename)
% INPUT
%  filename :  "boundary" file from Open Foam
% 
% OUTPUT
%  fp : n x 3 matrix of vital statistics for each physical boundary region
%     column 1: starting face
%     column 2: # of faces in each region
%     column 3: wall flag [1=wall, 0=not wall]
%  stats : original text lines from "boundary"
%     columns are: Name/nFaces/startFaces/WallFlag
%
fid = fopen(filename);

A = regexp(fileread(filename),'\n','split');
linenum = find(contains(A,'startFace'));
linenum2 = find(contains(A,'nFaces'));

% Will keep the NAME of the region
linenum3 = find(contains(A,'{'));
linenum3(1) = [];

% Finally, decide which regions are WALLS
linenum4 = find(contains(A,'type'));

%temp=textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum(1)-1)
%temp{1}

% I think this can be simplified. It should not
%   be necessary to reread the file, since we have each line stored in A{}
%   Come back to this later.
for i=1:length(linenum)
    f1(i) = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum(i)-1);
    frewind(fid);
    f2(i) = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum2(i)-1);
    frewind(fid);
    f3(i) = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum3(i)-2);
    frewind(fid);
    f4(i) = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum4(i)-1);
    frewind(fid);
end

% Maintain these lines of the file as strings
fstats = [f3.' f1.' f2.' f4.'];

% startFace, nFace
fp = [f1.' f2.'];
fp = cell2table(fp);

% Save these values as integers
fp.fp1 = erase(fp.fp1,["startFace", ";"]);
fp.fp2 = erase(fp.fp2,["nFaces", ";"]);

fp4=zeros(length(linenum4),1);
for j=1:length(linenum4)
    fp4(j)=strcmp(strtrim(erase(f4{j},["type",";"])),"wall");
end

fp.fp1 = str2double(fp.fp1);
fp.fp2 = str2double(fp.fp2);

% Add on "wall" flag
fp = table2array(fp);
fp = [fp fp4];

end