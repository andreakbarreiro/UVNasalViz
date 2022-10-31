function [centroids] = compCentroidsUV(r,face)
% 
% INPUTS:
%     r:            (nP x 2) UV-coordinates for each point
%     face:         (nF x 3) List of points associated w/ each face        

% OUTPUTS:  
%     centroids:    (nF x 2) centroids for each face

centroids = zeros(length(face),2);

for j1=1:length(face)
    
    temp            = r(face(j1,:),:);
    centroids(j1,:) = mean(temp);
  
end

% No need for any other processing; this is it!

end

