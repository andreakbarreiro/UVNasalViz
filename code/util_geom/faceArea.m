function s = faceArea(x,y,z)
% 
% xPts = [bdyPts(:,1) bdyPts(:,4) bdyPts(:,7)];
% yPts = [bdyPts(:,2) bdyPts(:,5) bdyPts(:,8)];
% zPts = [bdyPts(:,3) bdyPts(:,6) bdyPts(:,9)];

a = sqrt((x(:,1)-x(:,2)).^2+(y(:,1)-y(:,2)).^2+(z(:,1)-z(:,2)).^2);
b = sqrt((x(:,2)-x(:,3)).^2+(y(:,2)-y(:,3)).^2+(z(:,2)-z(:,3)).^2);
c = sqrt((x(:,1)-x(:,3)).^2+(y(:,1)-y(:,3)).^2+(z(:,1)-z(:,3)).^2);

p = (a+b+c)./2;

s = sqrt(p.*(p-a).*(p-b).*(p-c));
