function [l_b X_b Y_b Z_b]=findb(R,FACE)
% findb : find boundary of a set of faces
% INPUTS:
%    R      : (np x 3) (x,y,z) coordinates
%    FACE   : (nf x 3) List of points for each face
%
% OUTPUTS:
%    l_b    : List of line segments given in terms of indices for each
%    endpoint
%    X_b, Y_b, Z_b: x,y,z coordinates for those endpoints
%
%    

if size(R,2)==3
X=R(:,1);Y=R(:,2);Z=R(:,3);
else 
    X=R(:,1);Y=R(:,2);Z=X*0;
end
Num=length(FACE);l=zeros(3*Num,2);

for(i=1:Num)
    l(3*i-2:3*i,:)=[sort([FACE(i,1) FACE(i,2)]);
                    sort([FACE(i,1) FACE(i,3)]);
                    sort([FACE(i,2) FACE(i,3)])];
end
l=sortrows(l,[1 2]);

front=max(l(2:end-1,:)-l(1:end-2,:)~=0,[],2);
behind=max(l(2:end-1,:)-l(3:end,:)~=0,[],2);
combine=min([front behind],[],2);
mark=find(combine==1)+1;

if(max(l(1,:)~=l(2,:)))
    mark=[1;mark];
end
if(max(l(end,:)~=l(end-1,:)))
    mark=[mark;length(l)];
end
l_b=l(mark,:);
X_b=X(l_b)';Y_b=Y(l_b)';Z_b=Z(l_b)';

end