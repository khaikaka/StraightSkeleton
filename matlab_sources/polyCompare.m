function polyCompare(poly1,poly2,tol)
% Checks if two sets of polygons represent the same shape within Tol. Will
% fail for topological changes.
m=length(poly1);
m2=length(poly2);
if m~=m2
    error('polygons do not match')
elseif m==0
    return% matching empty polygons
end
for i=1:m
    poly1(i).xy=simplifyPoly(poly1(i).xy,tol);
    first1(i,:)=poly1(i).xy(1,:);
    poly2(i).xy=simplifyPoly(poly2(i).xy,tol);
    first2(i,:)=poly2(i).xy(1,:);
end
% order polygons by sorting starting points
[~,ind1]=sort(first1(:,1)-sqrt(2)*first1(:,2));
poly1=poly1(ind1);
[~,ind2]=sort(first2(:,1)-sqrt(2)*first2(:,2));
poly2=poly2(ind2);
for i=1:m
    if numel(poly1(i).xy)~=numel(poly2(i).xy) ||...
            ~all(abs(poly1(i).xy(:)-poly2(i).xy(:))<tol)
        error('polygons do not match')
    end
end

