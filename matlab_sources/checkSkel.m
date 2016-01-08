function halfedge=checkSkel(skel)

n=1:size(skel.vLRF,1);
halfedge=[n' skel.vLRF(:,1)
    n' skel.vLRF(:,2)
    n' skel.vLRF(:,3)];

checkManifold(halfedge);