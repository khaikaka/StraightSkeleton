function poly=closure(skel,h)
% ClOSURE - Polygon closure from a straight skeleton
%    Closure takes the vectorized skeleton, skel, of a polygon and uses
%    that to compute the outset of the inset of the original polygon by a
%    distance h. This results in a "rounded off" version of the original
%    polygon, similar to an erode/dilate from image processing. The output
%    polygon may contain a number of perimeters (wound CCW) and holes
%    (wound CW) in no particular order, having the form poly(i).xy.

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
        start=0;
        startLRF=0;
        last=skel.vLRF(current,LRF);
    elseif skel.z(next)<=h% new point
        if startLRF==0
            start=current;
            startLRF=LRF;
        end
        bR=skel.fCRL(current,LRF);
        bL=skel.fCRL(current,mod(LRF,3)+1);
        NR=skel.baseN(bR,:);
        NL=skel.baseN(bL,:);
        if NR*NL'<-.01 && NR(1)*NL(2)-NR(2)*NL(1)>0
            a=(skel.z(current)-h)/(skel.z(current)-skel.z(next));% convex combination
            xy=skel.xy(next,:)*a+skel.xy(current,:)*(1-a);
            NC=NR+NL;
            NC=NC/norm(NC);
            
            b=skel.xy(bR,:);
            c=skel.xy(skel.vLRF(bR,2),:);
            a=(NC-NR)*(xy-b)'/((NC-NR)*(c-b)');
            a=max(min(a,1),0);
            poly(p).xy(end+1,:)=c*a+b*(1-a);
            
            b=skel.xy(bL,:);
            c=skel.xy(skel.vLRF(bL,2),:);
            a=(NC-NL)*(xy-b)'/((NC-NL)*(c-b)');
            a=max(min(a,1),0);
            poly(p).xy(end+1,:)=c*a+b*(1-a);
        else
            a=skel.z(current)/(skel.z(current)-skel.z(next));% convex combination
            poly(p).xy(end+1,:)=skel.xy(next,:)*a+skel.xy(current,:)*(1-a);
        end        
        last=next;% switch direction to follow next face
    else% continue around face
        last=current;
        current=next;
    end
end