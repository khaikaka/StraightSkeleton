function polys=separate(loops,h,OUTSIDE)
% SEPARATE - organize polygons by containment
%    Separate takes a disorganized set of polygons, loops, of the form
%    loops(i).xy, where xy is an n-by-2 set of vertices. These polygons
%    should not intersect themselves or each other. The Boolean, OUTSIDE,
%    specifies if the polygons should all be inverted to represent the
%    negative space instead. If OUTSIDE is true, then h specifies the
%    margin the outer perimeter gets beyond the bounding box of the
%    polygons.
%
%    The output, polys, is of the form polys(i).poly(j).xy, where i is the
%    index of the polygon-with-holes, and j is the index of the polygon
%    within i. Any j==1 is a perimeter, wound CCW, while any j>1 is a hole,
%    wound CW. 
hpoly=sortPoly(loops);% create a heirarchical polygon
if OUTSIDE% expand bounding box for outsets
    lower=logical([1 1
        1 0
        0 0
        0 1]);
    hpoly.xy(lower)=hpoly.xy(lower)-3*h;
    hpoly.xy(~lower)=hpoly.xy(~lower)+3*h;
end
polys=flattenHPoly(hpoly,OUTSIDE);% split into separate polygons-with-holes
