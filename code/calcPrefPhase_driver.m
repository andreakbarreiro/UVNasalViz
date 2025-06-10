 % Driver to calculate phase preference



addpath('util_files');
addpath('util_geom');
addpath('util_plot');

% Inthavong directory
caseDir = '../caseDirs/intha_amp_50_T1p5/';

%% Read geometries, if not already done
read_polyMesh_flag = 0;
read_uv_flag       = 0;
if (read_polyMesh_flag)
    % I don't love having to keep multiple copies of the polyMesh
    % and UV files. However, it's probaly better than having endless
    % confusion
    %
    caseDirForPolymesh = '../caseDirs/intha_amp_50/';

    % Boundary file name
    fname   = 'constant/polymesh/boundary';

    [fp fstats]=facestats([caseDirForPolymesh fname]);
    [R FACE S lmax] = read_polyMesh(caseDirForPolymesh);
end

if (read_uv_flag)
    %% Load in boundary file
    bdydir      = '../../MeshCreation/ReadTecplot_GetRegions/'
    bdyfile     = 'bdy_All_Inthavong.mat'
    %bdyfile     = 'bdy_Olfac_Inthavong.mat';
    load([bdydir bdyfile]);

    %% This will contain the proper UV file name!!!
    [r face r2R]=read_uv(sourcefile);
    if (0)
        if (0)
            % Issue: UV coordinates are reversed, somehow, from uv2.obj
            % 
            % We appear to have the same number of faces, however. WHY???
            %
            uvfile = 'uv.obj';
            [r face r2R]=read_uv([caseDir uvfile]);
        end
        objdir  = '../../UVViz/caseDirs/Inthavong_OrthoOnly_0p5Hz/';
        objfile = 'uv2.obj';
    
        [r face r2R]=read_uv([objdir objfile]);
    end
end

%% Create phase map if, have not already done so
create_phase_map_flag=1;

%% Do I have warray in a file, or do I need to read it in
read_warray_from_file = 1;

if (create_phase_map_flag)
    if (read_warray_from_file)
        warrayFname = 'freq_0p64.mat';
        temp = load([caseDir warrayFname]);

        tarray = temp.tarray;
        nT     = length(tarray);
        warray = temp.warray;
        wallInfo = temp.wallInfo;;
        Tperiod = temp.Tperiod;
    else
        % Array of time points
        tarray = 0.01:0.01:5;  nT = length(tarray);
    
        % Set up storage for WSS
        nFace  = size(face,1);
        warray = zeros(nFace,nT);
    
        for j1=1:nT
            j1
            tstep = tarray(j1);
            wssfname    = sprintf('%s%g/wallShearStress',caseDir,tstep);
            [w,wallInfo] = read_wss_mag_OF_face(wssfname);
            warray(:,j1) = w;
        end
        %% Adjust for high flow rate during retronasal part of cycle?
        Tperiod     = 2.94;
    end

    adjustRetro = 0;
    if (adjustRetro)
        tarraySign = mod(tarray,Tperiod)-Tperiod/2;
        % For the retronasal parts of the cycle, divide WSS by this factor
        retroAdjustFactor = 1.68;
        adjustTimes = find(tarraySign*tarraySign(1)<0);

        warray(:,adjustTimes)= warray(:,adjustTimes)/retroAdjustFactor;
    end

    % Average over multiple periods
    [avgW,timeMod]      = calcAvgOverCycles(warray,tarray,Tperiod);

    %% Now try smoothing in space
    smoAvgW = mapQuantFaceToPoint(avgW,face,wallInfo);
    % Map back to faces
    faceSmoAvgW = mapQuantPointToFace(smoAvgW,face,wallInfo);

    % Get preferred phase for each point in space
    [prefPhase,maxOT]   = calcPrefPhase(avgW,timeMod);

    %% Plot preferred phase map
    % Get preferred phase for each point in space
    % For SMOOTHED data
    [prefPhaseSmo,maxOTS]   = calcPrefPhase(faceSmoAvgW,timeMod);

    %% Save the file, along with other parameters used
    if (adjustRetro)
        outfname = [caseDir 'prefPhaseSmo_retroAdjust.mat'];
    else
        outfname = [caseDir 'prefPhaseSmo_noAdjust.mat'];
    end
    save(outfname,'prefPhaseSmo',"maxOTS");
    save(outfname,'Tperiod', 'timeMod','tarray','-append');
    if (adjustRetro)
        save(outfname,'retroAdjustFactor','-append');
    end
end
