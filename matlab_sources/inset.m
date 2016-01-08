function poly=inset(skel,h)
% INSET - Polygon inset from a straight skeleton
%    Inset takes the vectorized skeleton, skel, of a polygon and uses that
%    to compute the polygon inset of the original polygon. The inset
%    distance is specified as h >= 0. The output polygon may contain a
%    number of perimeters (wound CCW) and holes (wound CW) in no particular
%    order, having the form poly(i).xy.

active=skel.z>h;
start=0;
startLRF=0;
current=find(active,1);
if isempty(current)
    poly=[];
    return
end
LRF=2;
last=skel.vLRF(current,LRF);
active=repmat(active,1,3);
active(current,LRF)=false;
p=1;
poly(p).xy=zeros(0,2);
poly(p).Vup=[];
poly(p).LRF=[];
while 1
    LRF=find(skel.vLRF(current,:)==last,1);
    next=skel.vLRF(current,mod(LRF,3)+1);
    active(current,LRF)=false;
    if current==start && LRF==startLRF% closed curve
        [current,LRF]=find(active,1);
        if isempty(current)% finished
            break
        end
        % start new polygon
        p=p+1;
        poly(p).xy=[];
        poly(p).Vup=[];
        poly(p).LRF=[];
        start=0;
        startLRF=0;
        last=skel.vLRF(current,LRF);
    elseif skel.z(next)<=h% new point
        if startLRF==0
            start=current;
            startLRF=LRF;
        end
        a=(skel.z(current)-h)/(skel.z(current)-skel.z(next));% convex combination
        a=max(min(a,1),0);% hedge against rounding error
        poly(p).xy(end+1,:)=skel.xy(next,:)*a+skel.xy(current,:)*(1-a);
        poly(p).Vup(end+1,:)=current;
        poly(p).LRF(end+1,:)=mod(LRF,3)+1;
        last=next;% switch direction to follow next face
    else% continue around face
        last=current;
        current=next;
    end
end