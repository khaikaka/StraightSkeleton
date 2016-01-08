function Tadj=adjacency(f)

nf=size(f,1);
halfedge=[f(:,[1 2])
    f(:,[2 3])
    f(:,[3 1])];% list of edges
ind1=checkManifold(halfedge);
ne=length(ind1);

Tadj=zeros(size(f));% triangle adjacency
for i=1:ne/2
    t1=mod(ind1(2*i-1)-1,nf)+1;
    e1=ceil(ind1(2*i-1)*3/ne);
    t2=mod(ind1(2*i)-1,nf)+1;
    e2=ceil(ind1(2*i)*3/ne);
    Tadj(t1,e1)=t2;
    Tadj(t2,e2)=t1;
end
