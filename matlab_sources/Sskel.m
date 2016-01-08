function skel=Sskel(poly,cutoff,DEBUG)
% Sskel - The straight skeleton of the polygon, poly.
%    The straight skeleton can be thought of as the dual of all the polygon
%    insets, and allows those insets to be rapidly calculated. It can also
%    be thought of as the roof shape, if the polygon were the walls of a
%    building, to make all the roof slopes identical. The Z-value of the
%    straight skeleton vertices mark the local half-width of the polygon,
%    hence the straight skeleton is useful for determining thickness.
%
%    The input polygon, poly, must be in the form poly(i).xy where the
%    first column is X-values and the second column is Y-values. poly(1)
%    must be the perimeter, wound CCW, while poly(k) for k > 1 must all be
%    holes contained in the perimeter, wound CW. Multiple perimeters will
%    result in an error. 
%
%    Polygons should not self-intersect. Mostly this will result in a
%    self-intersection error, but it is also possible to get a complete
%    straight skeleton, though in this case the insets may well also be
%    self-intersecting.
%
%    Sharp reflex angles are clipped when they exceed the specified angle,
%    cutoff, in degrees. Cutoff is limited to between 90 and 179 degrees,
%    with a default of 91.
%
%    Set the DEBUG flag to true to see a step-by-step plot of the
%    algorithm's progress. This is a good mode to set breakpoints in to
%    inspect for errors.
%
%    The output, skel, is a vectorized skeleton. skel.xy is m-by-2, where m
%    is the number of skeleton vertices, including those of the input
%    polygons. skel.z is m-by-1, skel.vLRF is m-by-3, containing the
%    indicies of the three vertices each one is attached to, in CCW order.
%    skel.fLRF points to the the left base-vertex of each of the three
%    faces that each vertex connects to. skel.fLRF(i,1)==0 for any i that
%    is base-vertex, since these points only attach to two faces.
%    skel.baseRefPoly is n-by-2, where n is the number of base-vertices,
%    referring to the first n skeleton points. skel.baseRefPoly(i,:)=[a,b],
%    meaning skel.xy(i,:)=poly(a).xy(b,:), indicating which input points
%    were kept. skel.baseN is also n-by-2, giving the unit-normal vectors
%    of each edge to the right of the corresponding base-vertex (pointing
%    away from the polygon). 

global bb scale
warning('OFF','all');% ignore matlab singular matrix warnings

if nargin<3
    DEBUG=false;
end
if nargin<2
    cutoff=91;
else
    cutoff=max(min(cutoff,179),90);
end

edgeX=[];
edgeY=[];
bottomV=[];
ref=[];
k=0;% polygon hole index (not including degenerates)
for j=1:length(poly)% process each polygon
    polyj=poly(j).xy;
    [polyj,refj]=deSliver(polyj,100*eps);% remove degenerate points
    n=size(polyj,1);
    if n==0% ignore degenerate polygons
        continue
    end
    ds=diff(polyj([1:end 1],:));
    N=ds*[0 -1;1 0];% normal vectors
    N=bsxfun(@rdivide,N,sqrt(sum(N.^2,2)));% normalize
    convex=N(:,1).*N([end 1:end-1],2)-N(:,2).*N([end 1:end-1],1)<0;
    flat=sum(N.*N([end 1:end-1],:),2)>cosd(cutoff);
    sharp=~convex & ~flat;% cut off these corners
    ind=1:n;
    ind1=ind(sharp);
    NR=N(ind1,:);
    ind=[ind(end) ind(1:end-1)];
    NL=N(ind(sharp),:);
    sharp=ind1+(0:length(ind1)-1);
    for i=1:length(sharp)% add edges to sharp reflex corners
        polyj=polyj([1:sharp(i) sharp(i) sharp(i)+1:end],:);
        refj=refj([1:sharp(i) sharp(i) sharp(i)+1:end]);
        Ni=NL(i,:)+NR(i,:);
        N=[N(1:sharp(i)-1,:)
            Ni/norm(Ni)
            N(sharp(i):end,:)];
    end    
    n=size(polyj,1);
    % create initial vertex objects
    L=size(edgeX,1);
    for i=1:n
        V(i+L)=vertex([],[],polyj(i,:));
    end
    % finish linking
    for i=1:n
        linkV(V(i+L),V(mod(i,n)+1+L));
        F=face(V(i+L),N(i,:));
        linkVF(V(i+L),F);
    end
    % calculate bisectors
    for i=1:n
        bisectBaseV(V(i+L));
    end
    if k>0% skip first polygon (outer) and find bottoms of holes
        [~,bottomV(k)]=min(polyj(:,2));
        bottomV(k)=bottomV(j-1)+L;
    end
    % save edge lists to calculate keyholes
    edgeX=[edgeX;polyj(:,1) polyj([2:end 1],1)];
    edgeY=[edgeY;polyj(:,2) polyj([2:end 1],2)]; 
    ref=[ref
        j*ones(size(refj))' refj'];
    k=k+1;
end
if isempty(edgeX)% return empty skeleton if polygon degenerates
    skel.xy=zeros(0,2);
    skel.z=zeros(0,1);
    skel.vLRF=zeros(0,3);
    skel.fLRF=zeros(0,3);
    return
end
n=length(V);
% global bounding box (for plotting) and scale (to calculate uncertainty)
bb=[min(edgeX(:,1)),max(edgeX(:,1)),min(edgeY(:,1)),max(edgeY(:,1))];
scale=max(abs(bb));
if DEBUG
    plot(V,'b');
    showNeighbors(V(1));
end
% calculate keyhole locations
for j=1:length(bottomV)
    % find edges that intersect a vertical line from bottomV(j)
    x=edgeX(bottomV(j),1);
    y=edgeY(bottomV(j),1);
    ind=find(edgeX(:,1)<=x & edgeX(:,2)>=x);
    a=(x-edgeX(ind,1))./(edgeX(ind,2)-edgeX(ind,1));
    p=[edgeX(ind,2).*a+edgeX(ind,1).*(1-a),...
        edgeY(ind,2).*a+edgeY(ind,1).*(1-a)];% convex combination
    [py,i]=sort(p(:,2),'descend');
    k=i(find(py<y,1));% find the nearst intersection below
    if isempty(k)
        error('hole outside perimeter')
    end
    lowerV=V(ind(k));% the left vertex of the intersected edge
    % insert keyholes
    Vnew=keyhole(lowerV,V(bottomV(j)),p(k,:));
    V=[V Vnew];% add new vertices
    added(j)=n+3*(j-1);% and keep track of where they are
    % update edge list with new vertices for next keyhole
    RV=[Vnew.RightV];
    LVxy=vertcat(Vnew.xy);
    RVxy=vertcat(RV.xy);
    edgeX(end+(1:3),:)=[LVxy(:,1) RVxy(:,1)];
    edgeY(end+(1:3),:)=[LVxy(:,2) RVxy(:,2)];
    edgeX(ind(k),2)=p(k,1);
    edgeY(ind(k),2)=p(k,2);
    if DEBUG
        plot(V,'b');
        showNeighbors(V(1));
        xy=vertcat(V.xy);
        uv=vertcat(V.bi);
        quiver(xy(:,1),xy(:,2),uv(:,1),uv(:,2),'g.')
%         keyboard
    end
end
botV=V(bottomV);% original keyhole vertices (not added)

% initialize the priority queue and begin the main loop
Q=vertex.empty;
for i=1:length(V);
    Q(end+1)=vertex(V(i),V(i).RightV);
    addQ(Q(end));
end
[~,q]=sort([Q.z]);% order priority queue
Q=Q(q);
while 1
    removeQ(Q(1));% remove this as a candidate from its parents
    % check that queue point is valid, with parents that are still
    % neighbors and no nearer candidates
    if ~isConnected(Q(1)) || ~isFirst(Q(1))
        Q(1)=[];% not a valid vertex        
        continue
    end
    
    % insert queue point into skeleton
    if ~isempty(Q(1).toRemove)
        for badV=Q(1).toRemove
            removeV(badV);% remove vertices beyond intersection
        end
    end    
    addV(Q(1));% attach new skeleton vertex
    V(end+1)=Q(1);% save new skeleton vertex
    Q(1)=[];% update priority queue
    
    % check for local termination
    if V(end).LeftNeighbor.LeftNeighbor==V(end)
        finalLink(V(end),V(end).LeftNeighbor);% connect missing link      
        if isempty(botV)% global termination
            if DEBUG
                plot(V,'b')
            end
            break
        else% open a keyhole and reset algorithm
            unKeyhole(botV(1),V(added(1)+(1:3)));  
            % reinitialize priority queue
            Q=vertex.empty;
            current=botV(1);
            i=1;
            while 1
                i=i+1;
                next=current.RightNeighbor;
                Q(end+1)=vertex(current,next);
                if Q(end).status==1
                    addQ(Q(end));
                end
                current=next;
                if current==botV(1)
                    break
                end
            end           
            [~,q]=sort([Q.z]);% order priority queue
            Q=Q(q);
            botV(1)=[];
            added(1)=[];
        end     
%         DEBUG=1;%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot progress, if desired
        if DEBUG
            plot(V,'b')
            showNeighbors(Q(1).LeftV)
            xy=vertcat(Q.xy);
            plot(xy(:,1),xy(:,2),'ko')
            V3=[Q.RightV];
            xy=vertcat(V3.xy);
            uv=vertcat(V3.bi);
            quiver(xy(:,1),xy(:,2),uv(:,1),uv(:,2),'g.')
%             keyboard
        end
        continue
    end
    
    % find new intersection on left
    i=1;
    clear VLeft
    VLeft(i)=V(end).LeftNeighbor;
    PVLeft=vertex(VLeft(i),V(end));% left possible vertex
    statusL=PVLeft.status;
    % continue around face to find valid intersection
    while statusL~=1 && ~VLeft(i).isbase% intersection is behind left bisector
        i=i+1;
        VLeft(i)=VLeft(i-1).RightV;
        permuteV(VLeft(i),VLeft(i-1));
        permuteV(VLeft(i-1).LeftV,VLeft(i-1));
        PVLeft=vertex(VLeft(i),V(end));
        statusL=PVLeft.status;
        if PVLeft.dL>VLeft(i-1).dR;
            statusL=2;% invalid
        end
    end
    
    % find new intersection on right
    i=1;
    clear VRight
    VRight(i)=V(end).RightNeighbor;
    PVRight=vertex(V(end),VRight(i));% right possible vertex
    statusR=PVRight.status;
    % continue around face to find valid intersection
    while statusR~=1 && ~VRight(i).isbase% intersection is behind right bisector
        i=i+1;
        VRight(i)=VRight(i-1).LeftV;
        permuteV(VRight(i),VRight(i-1));
        permuteV(VRight(i-1).RightV,VRight(i-1));
        PVRight=vertex(V(end),VRight(i));
        statusR=PVRight.status;
        if PVRight.dR>VRight(i-1).dL;
            statusR=3;% invalid
        end
    end
    
    % save valid intersections
    if statusL==2 && statusR==3
        error('self intersection')
    end  
    if statusL==1
        addQ(PVLeft);
        flagRemoval(PVLeft,VLeft(1:end-1));
        Q=sortInto(Q,PVLeft);
    end
    if statusR==1
        addQ(PVRight);
        flagRemoval(PVRight,VRight(1:end-1));
        Q=sortInto(Q,PVRight);
    end      
    
    % plot progress, if desired
    if DEBUG
        plot(V,'b')
        showNeighbors(V(end))
        if statusL==1
            plot(PVLeft.xy(1),PVLeft.xy(2),'ko')
        end
        if statusR==1
            plot(PVRight.xy(1),PVRight.xy(2),'ko')
        end
        V3=[V(end)
            V(end).RightNeighbor
            V(end).LeftNeighbor];
        xy=vertcat(V3.xy);
        uv=vertcat(V3.bi);
        quiver(xy(:,1),xy(:,2),uv(:,1),uv(:,2),'g.')
%         keyboard
    end 
end

% assign indices now that the skeleton is static
j=1;
for i=1:length(V)
    if V(i).status==1% only keep valid vertices
        index(V(i),j);
        Vout(j)=V(i); 
        j=j+1;
    end
end

% build vector skeleton
skel.xy=vertcat(Vout.xy);
skel.z=vertcat(Vout.z);
for j=1:length(Vout)    
    skel.vLRF(j,:)=[Vout(j).LeftV.ind Vout(j).RightV.ind Vout(j).ForwardV.ind];
    if isempty(Vout(j).CollapsedF)
        skel.fCRL(j,:)=[0 Vout(j).RightF.V.ind Vout(j).LeftF.V.ind];
    else
        skel.fCRL(j,:)=[Vout(j).CollapsedF.V.ind ...
            Vout(j).RightF.V.ind ...
            Vout(j).LeftF.V.ind];
    end
end
skel.baseRefPoly=ref;
for j=1:size(ref,1)%basev
    skel.baseN(j,:)=Vout(j).RightF.N;
end

% check graph validity
checkSkel(skel);
warning('ON','all');% reset warnings to default

function Q=sortInto(Q,next)
% should implement as proper priority queue to reduce from O(n) to
% O(log(n))
for j=1:length(Q)
    if next.z<=Q(j).z
        Q=[Q(1:j-1) next Q(j:end)];
        return
    end
end
Q=[Q next];

