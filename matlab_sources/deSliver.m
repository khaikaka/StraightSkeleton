function [poly,ind]=deSliver(poly,tol,odd,last,ind)
% Recursive function to remove slivers from the input polygon, poly, an 
% n-by-2 matrix of x-y coordinates. Slivers are edges folded
% back on themselves within the input, tol). When calling this function,
% only the first two arguments are specified. The output polygon, poly, is
% m-by-2, where m < n and ind is an m-by-1 vector corresponding to the
% retained vertices of the input poly.
n=size(poly,1);
if nargin<3% initialize for first call
    odd=false;
    last=true;
    ind=1:n;
end
if n<3% two-point polygon degenerates to nothing
    poly=[];
    return
end
% indices for operating on every-other edge pair
p1=odd+1:2:n-1;
p2=p1+1;
p3=mod(p2,n)+1;
% find folded edge pairs
v1=poly(p3,:)-poly(p2,:);
v2=poly(p2,:)-poly(p1,:);
d=v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);% cross product
L=max(sqrt([sum(v1.^2,2) sum(v2.^2,2)]),[],2);
colinear=abs(d)<=tol*L;
folded=sum(v1.*v2,2)<=tol*L;
% remove folded edge pairs
if any(colinear & folded)
    poly(p2(colinear & folded),:)=[];
    ind(p2(colinear & folded))=[];
    next=true;
else
    next=false;
end
% recursively continue until nothing has been removed on the last two times
if last || any(colinear & folded)
    [poly,ind]=deSliver(poly,tol,~odd,next,ind);
end