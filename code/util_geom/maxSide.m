function l = maxSide(x,y,z)

a = sqrt((x(:,1)-x(:,2)).^2+(y(:,1)-y(:,2)).^2+(z(:,1)-z(:,2)).^2);
b = sqrt((x(:,2)-x(:,3)).^2+(y(:,2)-y(:,3)).^2+(z(:,2)-z(:,3)).^2);
c = sqrt((x(:,1)-x(:,3)).^2+(y(:,1)-y(:,3)).^2+(z(:,1)-z(:,3)).^2);

l = max(max([a,b,c]));

end