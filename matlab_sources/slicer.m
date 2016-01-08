function s=slicer(model,z,wind)

v=model.v;
f=model.f;
if isfield(model,'Tadj')
    Tadj=model.Tadj;
elseif ~isempty(f)
    Tadj=adjacency(f,v);
end

ind=1:3;
s=struct('loop',cell(length(z),1));% slice structure
for i=1:length(z)% slice index
    vz=sign(v(:,3)-z(i));% above or below slice?
    vz(vz==0)=1;% exactly on slice is considered above slice
    tz=abs(sum(vz(f),2))<3;% triangles crossing this z-slice
    if nnz(tz)==0% empty slice
        s(i).loop=struct('xy',cell(1,0),'wind',cell(1,0));
        continue
    end
    loopi=1;% loop index
    
    tri=find(tz,1);% first triangle
    z1=vz(f(tri,:));
    edge=[ind(diff([z1;z1(1)])<0) ind(diff([z1;z1(1)])>0)];% ccw
    tri=[tri tri Tadj(tri,edge(end))];% next triangles
    tz(tri)=0;% visited
    while 1
        elast=ind(Tadj(tri(end),:)==tri(end-1));% last edge on current triangle
        enext=mod(elast+vz(f(tri(end),mod(elast-2,3)+1))-1,3)+1;% next edge (-vz for cw)
        tnext=Tadj(tri(end),enext);% next triangle
        if ~tz(tnext)% finished loop
            tri(end)=[];% remove last triangle (no edge)
            v1=v(f(sub2ind(size(f),tri,edge)),:);
            v2=v(f(sub2ind(size(f),tri,mod(edge,3)+1)),:);
            alpha=(z(i)-v1(:,3))./(v2(:,3)-v1(:,3));% convex parameter            
            alpha=min(max(alpha,0),1);% clamp ends to fix rounding error
            s(i).loop(loopi).xy=bsxfun(@times,1-alpha,v1(:,1:2))+bsxfun(@times,alpha,v2(:,1:2));
%             s(i).loop(loopi).wind=wind;% winding number
            % start next loop
            loopi=loopi+1;
            tri=find(tz,1);% first triangle
            if isempty(tri)% no more triangles in slice
                break
            end
            z1=vz(f(tri,:));
            edge=[ind(diff([z1;z1(1)])<0) ind(diff([z1;z1(1)])>0)];% ccw
            tri=[tri tri Tadj(tri,edge(end))];% next triangles
            tz(tri)=0;% visited
            elast=ind(Tadj(tri(end),:)==tri(end-1));
            enext=mod(elast+vz(f(tri(end),mod(elast-2,3)+1))-1,3)+1;% -vz for cw
            tnext=Tadj(tri(end),enext);
        end
        edge(end+1)=enext;
        tri(end+1)=tnext;
        tz(tnext)=0;% visited
    end
end