
% 
% UPDATE: After investigation
%   uv.obj is a perfect translation of the OpenFoam points
%      (modulo the "extra copy" of the boundary points)
%
%   uv2.obj has a different translation for wall vs. patch 
%      regions (nostrils, nasopharynx)
% The numerical values of the translations do not appear to be related to
% one another.
%
%
% QUESTION: HOW DO I FIX THIS WITH THE LEAST AMOUNT OF WORK
%   POSSIBLE?????
%  Here is what I will do.
%  - Write out new OBJ file.
%  -  AFAIK both uv and uv2 contain valid unwrappings, 
%       albeit with REVERSED ORIENTATIONS (i.e. left-right, I think, have
%       been swapped).
%  - FOR BOTH: remove extra "R" points
%  - FOR BOTH: find translation relative to OpenFoam points 
%  - FOR BOTH: wrote out "v" points ONLY
%     Since "vt" and "f" lines are unchanged: write these out
%     into a temporary file. Concatenate at the end!
%  - Concatenate in the file browser (i.e. outside of Matlab)
%   
addpath('../../UVViz/code/util_files');
addpath('../../UVViz/code/util_geom');
addpath('../../UVViz/code/util_plot');

% Let us try to revisit these few files
fdirlist = {'../../UVViz/caseDirs/Inthavong/';...
    '../../UVViz/caseDirs/Inthavong/';...
    '../../UVViz/caseDirs/intha_amp_50/';...
    '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/';...    
    };
fnamelist = {'Geometry_mapped_Orig.obj';...
    'Geometry_mapped_d4_Closed.obj';...
    'uv.obj';...
    'uv2.obj';...
    };

whF=4;
filename = [fdirlist{whF} fnamelist{whF}];

outfname = [fdirlist{whF} fnamelist{whF}(1:end-4) '_FIXED.obj']

fid=fopen(filename);

% We will open these temporary files, into which we write out
% the relevant lines
fid_v = fopen('out_v.txt','w');
fid_REST = fopen('out_REST.txt','w');

while ~feof(fid)
    tline=fgetl(fid);
    if length(tline)>2
        if min(tline(1:2)=='v ')
            % Why was this "end-2", rather than end?
            fprintf(fid_v,'%s\n',tline(3:end));
        else
            % WRITE OUT ENTIRE LINE
            %fprintf(fid_vt,'%s\n',tline(4:end));
            fprintf(fid_REST,'%s\n',tline);
        end
    end
end
fclose(fid);fclose(fid_REST);
fclose(fid_v);

% These each contain single array; we can use the "load" command
%r=load('out_vt.txt');map=load('out_f.txt');

[r,face,r2R,R]=read_3D_from_OBJ([fdirlist{whF} fnamelist{whF}]);

% Get rid of "extra points"
maxRind = max(r2R);
R(maxRind+1:end,:)=[];


%R = load('out_v.txt');

% Now OpenFoam version
%% Read OF geometry, if not already done
read_polyMesh_flag = 1;
%read_uv_flag       = 0;
if (read_polyMesh_flag)
    % I don't love having to keep multiple copies of the polyMesh
    % and UV files. However, it's probaly better than having endless
    % confusion
    %
    % Inthavong directory
    caseDir = '../../UVViz/caseDirs/intha_amp_50/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp,fstats]=facestats([caseDir fname]);
    [R_of,FACE_of,S_of,lmax_of] = read_polyMesh(caseDir);
end

temp = R_of - R;
temp = round(temp,4);

%% Make sure this is 3 or 6 values
unique(temp)

%% IF SO: we can just write our OpenFoam 
% values. 
if (length(unique(temp))<10)
    fid=fopen(outfname,'w');
    fprintf(fid,'# FIXED UV file\n');
    fprintf(fid,'# Original file: %s\n',filename);

    for j1=1:maxRind
        % Write out point
        fprintf(fid,'v %f %f %f\n',R_of(j1,1),R_of(j1,2),R_of(j1,3));
    end
    fclose(fid);
else
    error('Problem: points do not appear to be a translation!!');
end

