function polys=flattenHPoly(hpoly,INVERT)
% Take a heirarchical polygon and outputs a list, polys, of
% polygons-with-holes, poly. The first element of poly is the perimeter,
% the following elements are holes. If the INVERT flag is set to true, the
% function returns the negative space of the polygons, starting with the
% bounding box as a perimeter and the original periemters as holes within
% it.

if nargin<2
    INVERT=false;
end
if isempty(hpoly)
    polys=[];
    return
end

if INVERT
    polys=nextFlat(hpoly);
    polys=invert(polys);
else
    polys=[];
    for i=1:length(hpoly.hpoly)
        polys=[polys nextFlat(hpoly.hpoly(i))];
    end
end

function polys=nextFlat(hpoly)
polys(1).poly(1)=hpoly;% perimeter
if isfield(hpoly,'hpoly')
    for i=1:length(hpoly.hpoly)
        polys(1).poly(i+1)=hpoly.hpoly(i);% holes
        if isfield(hpoly.hpoly(i),'hpoly')
            for j=1:length(hpoly.hpoly(i).hpoly)
                polys=[polys nextFlat(hpoly.hpoly(i).hpoly(j))];
            end
        end
    end
end

function polys=invert(polys)
for i=1:length(polys)
    for j=1:length(polys(i).poly)
        polys(i).poly(j).xy=flipud(polys(i).poly(j).xy);
    end
end