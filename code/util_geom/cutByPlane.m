function [visFACE,visface]=cutByPlane(R,FACE,r,r2R,plane)
X=R(:,1);Y=R(:,2);Z=R(:,3);x=r(:,1);y=r(:,2);
a=plane(1);b=plane(2);c=plane(3);d=plane(4);
% cutByPlane:  Cut domain by plane
%  
% INPUTS:  
%   R, FACE:    coordinates, face vectors for 3D domain
%   r, face:    coordinates face vectors for corresponding UV domain
%   r2R:        mapping from r->R
%   plane:      plane specified as an R^4 vector w/
%               3D normal vector +offset
%               i.e. ax+by+cz = d
%
% OUTPUTS: 
%   visFACE, visface:  [0,1] vector indicating which faces are visible 


% Identify which side of the plane each point is on
judge=a*X+b*Y+c*Z-d;
i=find(judge>=0);
if(length(i)~=0)
    judge(i)=1;
end
i=find(judge<0);
if(length(i)~=0)
    judge(i)=-1;
end

% Points on each face
judgeFACE = sum(judge(FACE),2);
% Counts up 1 or -1 depending on which 
% side of face each point lies on

% Keep all faces which have at least 2 points on 
%   right sidde
visFACE=judgeFACE>0;
visface=[];


% From getplane
if (0)
l=zeros(3*length(FACE),2);
for(i=1:length(FACE))
    l(3*i-2,:)=[FACE(i,1) FACE(i,2)];
    l(3*i-1,:)=[FACE(i,2) FACE(i,3)];
    l(3*i,:)=[FACE(i,1) FACE(i,3)];
end
lj=judge(l);
lj=sum(lj,2);
% Line segment intersects plane
i=find(lj==0);
if(length(i)==0)
    Xj=[];Yj=[];Zj=[];xj=[];yj=[];
end
if(length(i)~=0)
    lj=l(i,:); 
    
    for(i=1:length(lj))
        j=find(r2R==lj(i,1));xa=x(j);ya=y(j);
        j=find(r2R==lj(i,2));xb=x(j);yb=y(j);
        
        Xa=X(lj(i,1));Xb=X(lj(i,2));
        Ya=Y(lj(i,1));Yb=Y(lj(i,2));
        Za=Z(lj(i,1));Zb=Z(lj(i,2));
        X0=Xb-Xa;Y0=Yb-Ya;Z0=Zb-Za;
    
        % Find coordinates of intersection along line
        u=(d-a*Xa-b*Ya-c*Za)/(a*X0+b*Y0+c*Z0);
        Xj(i)=Xa+u*X0;Yj(i)=Ya+u*Y0;Zj(i)=Za+u*Z0;

        % Now for UV coordinates
        if(length(xa)==1 & length(xb)==1)
            v=(Xb-Xj(i))/(Xj(i)-Xa);
            k1=v/(1+v);k2=1/(1+v);
            xj(i)=k1*xa+k2*xb;yj(i)=k1*ya+k2*yb;
        end
    end
    Xj=Xj';Yj=Yj';Zj=Zj';xj=xj';yj=yj';
end

Rj=[Xj Yj Zj];rj=[xj yj];
end

end