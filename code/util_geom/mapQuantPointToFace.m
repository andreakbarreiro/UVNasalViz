function [qOnFace] = mapQuantPointToFace(qOnPoint,face,wallInfo)
% 
% INPUTS:
%     qOnPoint:  (nP x nT)  Some quantity evaluated on faces 
%     face:      (nF x 3) List of points associated w/ each face        
%     wallInfo: Do I need this?

qOnFace = zeros(length(face),size(qOnPoint,2));

for j1=1:length(face)
    % Get values at point indices     
    temp = qOnPoint(face(j1,:),:);
    if (j1==1)
        size(temp)
        size(mean(temp))
    end
    % This is identical to a linear interpolation 
    %  at centroid
    %
    qOnFace(j1,:) = mean(temp);
end

% No need for any other processing; this is it!

end

