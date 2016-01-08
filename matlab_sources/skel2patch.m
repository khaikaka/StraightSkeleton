function p=skel2patch(skel,h,PLOT)
%SKEL2PATCH   Create a patch object from 2D Striaght Skeleton
%   SKEL2STL('filename',skel) writes an STL file of the 'roof' created
%   by the straight skeleton of a polygon. The output object also contains
%   the original polygon as the bottom, which should result in manifold
%   output geometry. 
%
%   SKEL2STL(...,h) clips the roof at height h, still outputting a manifold
%   object.
%
%   SKEL2STL(...,PLOT) controls the plotting functionality. Default is
%   'true'.
%
%   Example:
%
%     [x,y]=hilbert(3);
%     y(1)=-0.5;
%     y(end)=-0.5;
%     poly.xy=flipud([x' y']);
%     skel=Sskel(poly,0,0);
%     skel2stl('hilbert.stl',skel,0.02);

if nargin<3
    PLOT=true;% default to plotting
end
if nargin<2
    h=inf;% default to complete roof
end

% calculate inset for top
top=inset(skel,h);
% map skeleton vertices to inset vertices
V=zeros(0,3);
VupLRF=zeros(size(skel.vLRF));
topLRF=VupLRF;
for i=1:length(top)
    for j=1:length(top(i).Vup)
        VupLRF(top(i).Vup(j),top(i).LRF(j))=j;
        topLRF(top(i).Vup(j),top(i).LRF(j))=i;
    end
end

% process each skeleton face
F=zeros(0,3);
n=size(skel.baseRefPoly,1);
for i=1:n
    start=i;
    current=start;
    next=skel.vLRF(start,3); 
    face=[skel.xy(start,:) skel.z(start)];
    done=false;
    while ~done% continue CCW around face
        V3=skel.vLRF(next,:);
        j=find(V3==current,1);
        above=VupLRF(next,j);
        abovei=topLRF(next,j);
        if above% higher than h
            face(end+1,:)=[top(abovei).xy(above,:) h];
            k=mod(above-2,length(top(abovei).Vup))+1;
            face(end+1,:)=[top(abovei).xy(k,:) h];
            current=top(abovei).Vup(k);
            next=skel.vLRF(current,top(abovei).LRF(k));            
        else
            face(end+1,:)=[skel.xy(next,:) skel.z(next)];
            current=next;            
            next=V3(mod(j+1,3)+1);            
        end       
        done=next==start;
    end
    xyz=face;  
    m=size(xyz,1);
    V1=zeros(0,3);
    for k=1:m% remove duplicate points
        if any(xyz(k,:)~=xyz(mod(k,m)+1,:))
            V1(end+1,:)=xyz(k,:);
        end
    end
    poly.xy=V1(:,1:2);
    [tri,V2]=triPoly(poly);
    for j=(size(V1,1)+1):size(V2,1)
        [~,k]=min(sum(bsxfun(@minus,V2(j,:),V1(:,1:2)).^2,2));
        tri(tri==j)=k;
    end
    F=[F
        tri+size(V,1)];
    V=[V
        V1];
end

% process each input polygon (bottom face)
for i=1:skel.baseRefPoly(end,1);
    xy=skel.xy(find(skel.baseRefPoly(:,1)==i),:);
    m=size(xy,1);
    V1=zeros(0,2);
    for k=1:m% remove duplicate points
        if any(xy(k,:)~=xy(mod(k,m)+1,:))
            V1(end+1,:)=xy(k,:);
        end
    end
    poly(i).xy=V1;
end
[tri,V2]=triPoly(poly);
F=[F
    tri(:,[2 1 3])+size(V,1)];% reverse normals
V=[V
    V2 zeros(size(V2,1),1)];

% process each inset polygon (top face)
[tri,V2]=triPoly(top);
F=[F
    tri+size(V,1)];
V=[V
    V2 h+zeros(size(V2,1),1)];

[V,~,j]=unique(V,'rows');% reindex for unique vertices
F=j(F);% reindex faces
try
    adjacency(F);% check manifold
catch
    warning('Nonmanifold output')
end

p.vertices=V;
p.faces=F;

if PLOT
    figure(1)
    patch(p,'facecolor','r')
    axis equal
    axis vis3d
    view(3)
end
