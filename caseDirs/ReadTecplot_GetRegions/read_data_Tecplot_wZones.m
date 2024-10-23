function [R,P,FACE,S,lmax,zo_Info]=read_data_Tecplot_wZones(filename)
% Read data from a Tecplot.360 file
% 
% See Tecplot 360 EX 2021 Release 2
% "Data Format Guide" for detailed description
% 
% _wZones: INCLUDE PHYSICAL ZONES
%

% Find out number of points of each part
fid=fopen(filename);

% Read entire file as a sequence of characters
scan=fscanf(fid,'%c');

% First, look for indicators that a list of coordinates is
%    to be printed
%
% pa_num = number of regions to be read in
i=findstr('N=',scan);pa_num=length(i);pt_num=zeros(1,pa_num);
k0=2;
% AKB
% In my version of 'DUN_...' it is 'Nodes=' not 'N='
if isempty(i)
    i=findstr('Nodes=',scan);pa_num=length(i);pt_num=zeros(1,pa_num);
    k0=6;
end

% Similarly, see how many times we named a zone
zoStarts = findstr('ZONE T=',scan);zo_num=length(zoStarts);
zo_Names = cell(zo_num,1);
z0=8;

if (zo_num ~= pa_num)
    % These should be the same
    warning('Numbers of zones is inconsistent');
end

% After this, pt_num(j) should contain the number of points in
%    region j
for j=1:pa_num
    k=k0;
    while scan(i(j)+k)~=','
        k=k+1;
    end
    pt_num(j)=str2num(scan(i(j)+k0:i(j)+k));
end

% After this, zo_Names should contain the zone names
for j=1:pa_num
    k=z0;
    while (scan(zoStarts(j)+k)~= ' ')
        k=k+1;
    end
    zo_Names{j} = scan(zoStarts(j)+z0:zoStarts(j)+k-1);
end

%zo_Names

zo_Info = {zo_Names,pt_num'};

fclose(fid);

fidin=fopen(filename);
fidout=fopen('out.txt','w');
while ~feof(fidin)
    % Trim off leading whitespace
    tline=strtrim(fgetl(fidin));
    if length(tline)>0
        if ~isletter(tline(1))
            fprintf(fidout,'%s\n',tline);
        end
    end
end
fclose(fidin);fclose(fidout);

data1=importdata('out.txt');

% AKB: behavior of importdata has changed since this 
%    GUI was written??
% Actual matrix is part of a structure
data=data1.data;

%delete out.txt;

stat=sum(~isnan(data),2);
stat1=stat(1:end-1);stat2=stat(2:end);
res=(stat2+stat1)+(stat2-stat1)/10;
cut=find(res==8.2);  %Goes from 3->5
% AKB: This is identifying the last line which lists a face
% (which here is always a triple)
for i=1:length(cut)
    a=data(cut(i),1);
    % This is a check to ensure what we are working with is an INTEGER
    a=a-fix(a);   
    if a~=0
        cut(i)=nan;
    end
end
i=find(isnan(cut)==0);cut=cut(i);
cut=[0;cut;length(stat)];

plus=zeros(1,pa_num-1);
for(i=1:pa_num-1)
    plus(i)=sum(pt_num(1:i));
end

% Region 1: Saving the data for 3 coordinates + 1 data value, blockstyle
Num=pt_num(1);n=ceil(Num/5);
X=data(1:n,:)';
Y=data(n+1:2*n,:)';
Z=data(2*n+1:3*n,:)';
P=data(3*n+1:4*n,:)';
X=X(:);X=X(1:Num);
Y=Y(:);Y=Y(1:Num);
Z=Z(:);Z=Z(1:Num);
P=P(:);P=P(1:Num);

% Region 1: saving faces (as triples of integers, each identifies a point)
FACE=data(4*n+1:cut(2),1:3);

% Repeat, for regions 2-pa_num
if(pa_num>1)
    for(i=2:pa_num)
        Num=pt_num(i);n=ceil(Num/5);
        X1=data(cut(i)+1:cut(i)+n,:)';
        Y1=data(cut(i)+n+1:cut(i)+2*n,:)';
        Z1=data(cut(i)+2*n+1:cut(i)+3*n,:)';
        P1=data(cut(i)+3*n+1:cut(i)+4*n,:)';
 
        X1=X1(:);X1=X1(1:Num);
        Y1=Y1(:);Y1=Y1(1:Num);
        Z1=Z1(:);Z1=Z1(1:Num);
        P1=P1(:);P1=P1(1:Num);
        
        FACE1=data(cut(i)+4*n+1:cut(i+1),1:3)+plus(i-1);
        X=[X;X1];Y=[Y;Y1];Z=[Z;Z1];P=[P;P1];FACE=[FACE;FACE1];
    end
end

R=[X Y Z];
XXX=X(FACE);YYY=Y(FACE);ZZZ=Z(FACE);
R1=[XXX(:,1) YYY(:,1) ZZZ(:,1)];
R2=[XXX(:,2) YYY(:,2) ZZZ(:,2)];
R3=[XXX(:,3) YYY(:,3) ZZZ(:,3)];
S=faceArea(R1,R2,R3);

Num=length(FACE);
l=zeros(3*length(FACE),2);
for i=1:Num
    l(3*i-2,:)=[FACE(i,1) FACE(i,2)];
    l(3*i-1,:)=[FACE(i,2) FACE(i,3)];
    l(3*i,:)=[FACE(i,1) FACE(i,3)];
end
X1=X(l);Y1=Y(l);Z1=Z(l);
l=[X1(:,1) Y1(:,1) Z1(:,1)]-[X1(:,2) Y1(:,2) Z1(:,2)];
l=sqrt(l(:,1).^2+l(:,2).^2+l(:,3).^2);
lmin=min(l);lmax=max(l);

end