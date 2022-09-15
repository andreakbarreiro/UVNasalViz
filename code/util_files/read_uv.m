function [r face r2R]=read_uv(filename)
% read_uv: read "texture" vertices from an .obj file
%   Caution: when unwrapping, use "Keep original point order" option
%
%   We ASSUME this is a tetrahedral mesh: all faces must be
%     triangles
%  INPUTS:
%     filename:  Contains 3D model+UV-unwrapped model
%    
%
%  OUTPUTS:
%     r:        (np x 2 ) point coordinates (u,v)
%     face:     (nf x 3) point indices of each face
%     r2R:      (np x 1) map BACK to original 3D points R
%
% fid=fopen('model_unwrapped8.obj');
fid=fopen(filename);

% We will open these temporary files, into which we write out
% the relevant lines
fid_vt=fopen('out_vt.txt','w');fid_f=fopen('out_f.txt','w');
while ~feof(fid)
    tline=fgetl(fid);
    if length(tline)>2
        if min(tline(1:2)=='vt')
            fprintf(fid_vt,'%s\n',tline(4:end-2));
        end
        if min(tline(1:2)=='f ')
            i=find(tline=='/');
            tline(i)=' ';
            fprintf(fid_f,'%s\n',tline(3:end));
        end
    end
end
fclose(fid);fclose(fid_vt);fclose(fid_f);

% These each contain single array; we can use the "load" command
r=load('out_vt.txt');map=load('out_f.txt');
delete out_vt.txt; delete out_f.txt;

if (size(map,2)==9)
    % Face elements can be specified v/vt/vn or v/vt/vn
    %    if the latter, discard the "vn" as we don't use it.
    map(:,[3 6 9]) = [];    % AMS
end

% To map back to original 3D points
% Each row contains a pair of the form "orig point"/"new point"
%    indicating point IDs in the original "v" coordinates and the texture
%    "vt" coordinates respectively
%

r2R=[map(:,1:2);map(:,3:4);map(:,5:6)];

% We ASSUME there is no ambiguity; that each UV point maps to one and
%    only one XYZ point. However, an XYZ point may map to
%    more than one UV point. Therefore we sort on UV point index.
r2R=sortrows(unique(r2R,'rows'),2);

%tempSave = r2R;

r2R=r2R(:,1);

%p=P(r2R);

face=map(:,[2 4 6]);

end