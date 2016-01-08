function [tri,V]=triPoly(poly)
% Triangulates a single polygon with holes. Self intersection?

E=zeros(0,2);
V=zeros(0,2);
for i=1:length(poly)
    ind=(1:size(poly(i).xy,1))'+size(V,1);
    E(ind,:)=[ind ind([2:end 1])];% edges
    V(ind,:)=poly(i).xy;% vertices
end
if isempty(V)
    tri=[];
    return
end
DT=delaunayTriangulation(V,E);% constrained to polygon
tri=DT.ConnectivityList;
keep=isInterior(DT);
tri=tri(keep,:);
V=DT.Points;