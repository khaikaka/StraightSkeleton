function dxf2stl(filein,fileout,h,simpTol)

[~,c_Poly]=f_LectDxf(filein);
for i=1:size(c_Poly,1)
    loops(i).xy=simplifyPoly(flipud(c_Poly{i,1}),simpTol);
end
polys=separate(loops,h,0);

P.vertices=[];
P.faces=[];
for i=1:length(polys)
    skel=Sskel(polys(i).poly,179);
    p=skel2patch(skel,h,0);
    P.faces=[P.faces
        p.faces+size(P.vertices,1)];
    P.vertices=[P.vertices
        p.vertices];   
end

patch2stl(fileout,P)

figure(1)
clf
patch(P,'facecolor','r')
axis equal
axis vis3d
view(3)