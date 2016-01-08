function testSTLskeleton(filename)
close all

t=0.2;% layer thickness
simpTol=1e-2;
h=2;%exp(0.1);
VIZ=true;

model=loadSTL(filename);

phi=0;%pi/4;
% R=[cos(phi) sin(phi) 0
%     -sin(phi) cos(phi) 0
%     0 0 1];
R=[1 0 0
    0 cos(phi) sin(phi)
    0 -sin(phi) cos(phi)];
model.v=model.v*R;
model.v(:,3)=model.v(:,3)-min(model.v(:,3));% put on platform
model.n=model.n*R;

z=t/2:t:max(model.v(:,3))-t/2;% slice heights

tic
model.Tadj=adjacency(model.f);
toc

tic
smodel=slicer(model,z,1);
toc
% tic
% slices=sliceMesh(model.f,model.v,z);
% toc

tic
for j=1:length(smodel)
    j
    for i=1:length(smodel(j).loop)
        smodel(j).loop(i).xy=simplifyPoly(smodel(j).loop(i).xy,simpTol);
    end
    multiSkel(smodel(j).loop,h,VIZ);
end
toc

figure
% plotslice1(inter,z)
plotslice(smodel,z)

function []=multiSkel(loops,h,VIZ)
if VIZ
    figure(1)
    clf
    hold on
    axis equal
end
for OUTSIDE=[false,true];
    polys=separate(loops,h,OUTSIDE);
    for i=1:length(polys)
        skel=Sskel(polys(i).poly);
        if VIZ
            if OUTSIDE
                plotSkel(skel,'b')
            else
                plotSkel(skel,'r')
            end
        end
        in=inset(skel,h);
        if VIZ
            plotPoly(in,'k--')
        end
    end
end
drawnow