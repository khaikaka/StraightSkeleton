function a=testSkel()

tol=1e-10;
n=4;
a=[60
    94
    293.3005271015444
    100.2593587921374
    285.1946386414396
    35.114545799787429
    305.6865501127597
    53.745842001260691
    47.941562738578220
    196.9975907869691
    rand(10,1)*360
    ]';


% test vertex removal
poly1.xy=[
    0 1
    0 -1
    1 -1
    1 -3
    4 -3
    4 -1
    0.5 0
    4 1
    4 3
    1 3
    1 1];

skel=testRotation(poly1,a,tol,5);
b=inset(skel,1.5);
c=closure(skel,1.5);

% test manifold error checking
halfedge=checkSkel(skel);
halfedge(1,:)=halfedge(1,[2 1]);
try
    checkManifold(halfedge);
    error('Orientation check type 2 error')
catch
end
halfedge(1,1)=0;
try
    checkManifold(halfedge);
    error('Manifold check type 2 error')
catch
end
halfedge(1,:)=[];
try
    checkManifold(halfedge);
    error('Manifold check type 2 error')
catch
end

% test holes
poly1(2).xy=[2 1
    2 2
    3 2
    3 1];

poly1(3).xy=flipud([2 -1
    2 -2
    3 -2
    3 -1]);

% test debug mode
skel=Sskel(poly1,91,1);
plotSkel(skel,'b')
h=linspace(0,max(skel.z),4+2);
h=h(2:end-1);
% test empty polygons
in=inset(skel,2*h(end));
% test plotting
for i=1:length(h)
    in=inset(skel,h(i));
    if i==3
        plotPoly(in,'g--')
        clsr=closure(skel,h(i));
        plotPoly(clsr,'r')
    else
        plotPoly(in,'k--')
    end
end
drawnow

testRotation(poly1,a,tol,n);

% test polyCompare error checking
in2=in;
in2(end).xy(1,1)=in2(end).xy(1,1)+0.01;
try
    polyCompare(in,in2);
    error('polyCompare type 2 error')
catch
end
in2(end)=[];
try
    polyCompare(in,in2);
    error('polyCompare type 2 error')
catch
end

% test degenerate cases
poly1a=poly1(1);
poly1a.xy(7,:)=[];
testRotation(poly1a,a,tol,4);
poly1a.xy(1:2,1)=1;
testRotation(poly1a,a,tol,n);

poly5.xy=[0 0
    4 0
    3 1
    4 0
    4 2
    5 3
    5 4
    5 3
    4 2
    0 2
    0 3
    0 3
    0 3
    0 -1];
poly5(2).xy=[1 1
    2 1];
testRotation(poly5,a,tol,4);

poly6.xy=[-3 -3
    3 -3
    3 3
    -3 3];
poly6(2).xy=[1 1
    1 2
    2 2
    2 1];
poly6(3).xy=poly6(2).xy-3;
poly6(4).xy=poly6(2).xy;
poly6(4).xy(:,1)=poly6(4).xy(:,1)-3;
poly6(5).xy=poly6(2).xy;
poly6(5).xy(:,2)=poly6(5).xy(:,2)-3;
testRotation(poly6,[90 -90 180],tol,4);

poly7.xy=[0 0
    2.5 0
    2.5 3
    2 2.9
    1.8 3.5
    2 4.1
    2.5 4
    2.5 7
    0 7
    0 4
    0.5 4.1
    .7 3.5
    0.5 2.9
    0 3];
testRotation(poly7,a,tol,4);

poly8.xy=[-1 -1
    1 -1
    1 1
    1 -1];
Sskel(poly8);

% test self-intersection
poly1(2).xy(:,2)=poly1(2).xy(:,2)-1.5;
try
    Sskel(poly1);
    error('failed to find self intersection')
catch
end
poly1(3).xy(:,1)=poly1(3).xy(:,1)+2.5;
try
    Sskel(poly1);
    error('failed to find exterior hole')
catch
end

% test side-by-side holes
poly2.xy=[0 0
    1 1
    .1 2
    1 3
    .1 4
    2 7
    -2 7
    -.1 3
    -1 2
    -.1 1];

poly2(2).xy=[-0.05 5
    -1 6
    -0.05 6];

poly2(3).xy=[0.05 5
    0.05 6
    1 6];

skel2=testRotation(poly2,a,tol,n);

% test high resolution circle
theta=linspace(-pi,pi,100)';
theta=theta(2:end);
poly3.xy=[cos(theta) sin(theta)];

skel3=testRotation(poly3,.001,tol,n);

% test fractal symmetry
[x,y]=hilbert(3);
y(1)=-0.5;
y(end)=-0.5;
poly4.xy=flipud([x' y']);

skel4=testRotation(poly4,.001,tol,4);


function skel=testRotation(poly,r,tol,n)
skel=Sskel(poly);
h=linspace(0,max([skel.z]),n+2);
h=h(2:end-1);
for i=1:length(h)
    ins(i).in=inset(skel,h(i));
end
clsr=closure(skel,h(end));
for a=r
    R=[cosd(a) sind(a);-sind(a) cosd(a)];
    poly1=rotatePoly(poly,R);
    skel1=Sskel(poly1);   
    for i=1:length(h)            
        in1=inset(skel1,h(i));
        in1=rotatePoly(in1,R');
        polyCompare(ins(i).in,in1,tol);
    end
    clsr1=closure(skel,h(end));
    polyCompare(clsr,clsr1,tol);
end