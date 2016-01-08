function poly=simplifyPoly(poly,tol)
% Removes duplicate and colinear points within Tol. Input/output poly is a
% single polygon, n-by-2 matrix, where the two columns are X and Y values,
% respectively.
[~,start]=min(poly(:,1)-sqrt(2)*poly(:,2));% arbitrary order to avoid degeneracies
poly=circshift(poly,1-start);
i=1;
j=2;
k=3;
n=size(poly,1);
keep=true(n,1);
while k~=2
    if isColinear(poly(i,:),poly(j,:),poly(k,:),tol)
        keep(j)=false;     
    else
        i=j;
    end
    j=k; 
    k=mod(k,n)+1;
end
poly=poly(keep,:);

function colinear=isColinear(p1,p2,p3,tol)
v1=p3-p1;
v2=p2-p1;
d=v1(1)*v2(2)-v1(2)*v2(1);% cross product
colinear=abs(d)<tol*sqrt(sum(v1.^2));