function hpoly=sortPoly(poly)
% Takes an arbitrary set of geometrically valid polygons (winding number
% should only be zero or one for any point in R2) and sorts them into a
% heirarchy, hpoly, based on containment. Going one level down switches
% from perimeters (CCW) to holes (CW), or vice versa. 

if isempty(poly)
    hpoly=[];
    return
end
E=zeros(0,3);
V=zeros(0,2);
n=length(poly);
for i=1:n
    m=size(poly(i).xy,1);
    ind=(1:m)'+size(V,1);
    E(ind,:)=[ind ind([2:end 1]) i*ones(m,1)];% edges
    V(ind,:)=poly(i).xy;% vertices
    p(i,:)=poly(i).xy(1,:);% first point in polygon
end

inpoly=false(n);
for i=1:n
    R1=V(E(:,1),1)>p(i,1);
    R2=V(E(:,2),1)>p(i,1);
    E1=E(xor(R1,R2),:);% crossing x
    D1=V(E1(:,1),2)>p(i,2);
    D2=V(E1(:,2),2)>p(i,2);
    keep1=all([D1 D2],2);% above
    E2=E1(xor(D1,D2),:);% check for which side
    if ~isempty(E2)
        a=(V(E2(:,1),1)-p(i,1))./(V(E2(:,1),1)-V(E2(:,2),2));
        a=max(min(a,1),0);% clamp a to [0,1]
        keep2=(1-a).*V(E2(:,1),2)+a.*V(E2(:,2),2)>p(i,2);
    else
        E2=[];
        keep2=[];
    end
    E1=[E1(keep1,:);E2(keep2,:)];
    inpoly(i,:)=logical(mod(histc(E1(:,3),1:n),2))';
    inpoly(i,i)=0;% polygons are not inside themselves
end

minmax=[min(V);max(V)];
mmx=minmax(:,1);
mmy=minmax(:,2);
hpoly.xy=[mmx([1 1 2 2]) mmy([1 2 2 1])];% bounding box

hpoly.hpoly=nextLevel(poly,inpoly,1:n);


function hpoly=nextLevel(poly,inpoly,ind)

thisLevel=~any(inpoly,2);
if ~all(sum(inpoly(~thisLevel,ind(thisLevel)),2)==1)
    error('self intersection')
end
k=0;
poly(1).hpoly=[];
for i=ind(thisLevel)
    k=k+1;
    hpoly(k)=poly(i);
    under=~thisLevel & inpoly(:,i);
    inpoly1=inpoly(under,:);
    inpoly1(:,i)=0;
    if nnz(under)>0
        hpoly(k).hpoly=nextLevel(poly,inpoly1,ind(under));
    end
end
