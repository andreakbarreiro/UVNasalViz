function [temp] = mapQuantFaceToPoint(qOnFace,face,wallInfo)
% 
% INPUTS:
%     qOnFace:  (nF x nT) Some quantity evaluated on faces 
%     face:     (nF x 3) List of points associated w/ each face        
%     wallInfo: Do I need this?

% All points we need to be concerned about
allPoints = unique(reshape(face,[3*length(face),1]));


qOnPoints = zeros(size(allPoints,1),size(qOnFace,2));

size(qOnPoints)

for j1=1:length(allPoints)
    % Find all faces contain this point 
    [m n]=find(face == allPoints(j1));
%     if (j1==1)
%         size(qOnFace(m,:))
%         size(mean(qOnFace(m,:)))
%     end
    qOnPoints(j1,:) = mean(qOnFace(m,:));
       
end

% We now map back to original point indices 
%  (which may be identical given how we read in the data)
%  Use NaN for missing data
maxPointInd = max(allPoints);

temp   = nan(maxPointInd,size(qOnPoints,2));

% Each qOnPoint value should now be mapped to its correct index
temp(allPoints,:) = qOnPoints;


end

