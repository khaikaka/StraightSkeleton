function ind1=checkManifold(halfedge)

[halfedge,ind]=sort(halfedge,2);% make all edges go the same way
[~,~,j]=unique(halfedge,'rows');% find duplicate edges
[js,ind1]=sort(j);% sort the duplicates into pairs

% ensure all edges are paired
if mod(length(js),2)~=0
    error('Not Manifold (odd number of halfedges)')
end
check=reshape([diff(js);1],2,[])==0;
if ~all(check(1,:)) || ~all(~check(2,:))
    error('Not Manifold')
end
% ensure all pairs go opposite ways
check2=reshape(ind(ind1,1),2,[])==1;
if ~all(xor(check2(1,:),check2(2,:)))
    error('Not Oriented')
end