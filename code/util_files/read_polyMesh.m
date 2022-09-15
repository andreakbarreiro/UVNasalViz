function [p, sigFace, S, lmax] = read_polyMesh(caseDir)
% read_polyMesh: read points from an OpenFoam file
%   We ASSUME this is a tetrahedral mesh: all faces must be
%     triangles
%  INPUTS:
%     caseDir:    relative location of an OpenFoam case directory   
%        Must contain files: 
%           constant/polymesh/boundary
%           constant/polymesh/points      
%
%  OUTPUTS:
%     p:        (np x3) point coordinates (x,y,z)
%     sigFace:  (nf x3) point indices of each face
%     S:        (nfx1) area of each face
%     lmax:     unsure
%
    [fp fstats] = facestats([caseDir 'constant/polymesh/boundary']);

    % points from polymesh
    p = read_points([caseDir 'constant/polymesh/points']);
 
    % For each physical region
    startFaceIndex  =  fp(:,1) + 1;     %Offset by 1 because OpenFoam 
                                        %starts lists w/ "0"
    nFace           =  fp(:,2);
   
    % Read faces
    faces = read_faces([caseDir 'constant/polymesh/faces']);
    
    % Determine which boundary pieces are WALLS
    whichWalls = find(fp(:,3));
     
    % How many boundary faces do we have?
    nBdyFaces   = sum(nFace(whichWalls));    
    sigFace     = zeros(nBdyFaces,3);
    
    % Iterate through wall regions. 
    % Save ONLY boundary faces
    kk = [];
    firstI = 1;
    for j=1:length(whichWalls)        
        regI  = whichWalls(j);
        [regI nFace(regI)]
        % Position where to place coordinates in
        %   sigFace
        lastI  = firstI + nFace(regI)-1;
        
        % Indices for "faces"
        origFirstI  = startFaceIndex(regI);
        origLastI   = startFaceIndex(regI)+nFace(regI)-1;
        
        sigFace(firstI:lastI,:) = faces(origFirstI:origLastI,:);
        
        regI
        kk = [kk; regI firstI lastI origFirstI origLastI];
        
        % Reset for the next region
        firstI = firstI + nFace(regI);
        
    end
    %kk
    
    % Add "1" to the entire matrix, 
    % because OpenFoam starts indexing at "0"
    sigFace = sigFace + 1;
    
    % At this point we can clear faces
    clear faces
    
   % Actually, the below is NOT true
   % We do not need to interpolate to points
%     % Each row of sigFace contains three point indices
%     %   Our issue now is that we have plotting routines which 
%     %   require functions specified on POINTS, not necc. on FACES
%     %   So we need to be able to interpolate a FACE quantity to individual
%     %   points
%     % 
%     % This is a ham-handed way to do that:
%     %   We now get a list of WHICH points asoc. with each boundary face
%     %   
%     % On second thought, handle this later. 
%     %   For now cut down # of points
%     %
    sigFacePointsNoDup = unique(reshape(sigFace,[3*length(sigFace) 1]));
    length(sigFacePointsNoDup)
    sigFacePointsNoDup(1:10)
    sigFacePointsNoDup(end-10:end)
    
    % Discard any points AFTER this
    p = p(1:max(sigFacePointsNoDup),:);
    
    % Compute areas, lmax
    bdyPts = zeros(nBdyFaces,9);
    for j=1:nBdyFaces
        bdyPts(j,:) = reshape(p(sigFace(j,:),:).',[1 9]);
    end
    
    % x, y, z
    xPts = bdyPts(:,[1 4 7]);
    yPts = bdyPts(:,[2 5 8]);
    zPts = bdyPts(:,[3 6 9]); 
    
    % Area of each face
    S = faceArea(xPts,yPts,zPts);

    % max length of a side of a triMesh element
    lmax = maxSide(xPts,yPts,zPts);
end

