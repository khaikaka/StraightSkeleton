classdef vertex < handle
    
    properties(SetAccess = private)
        xy% horizontal position
        z=inf;% distance from associated edges (height of roof)
        delta% uncertainty, due to numerical precision
        LeftV=vertex.empty;% parent vertex (or neighboring base)
        RightV=vertex.empty;% parent vertex (or neighboring base)
        ForwardV=vertex.empty;% next vertex
        bi% bisector unit-vector towards ForwardV
        dL% distance ahead of left parent
        dR% distance ahead of right parent
        QL=vertex.empty;% list of ForwardV candidates on Left
        QR=vertex.empty;% list of ForwardV candidates on Right
        toRemove=vertex.empty;% vertices ahead of a parent
        LeftNeighbor=vertex.empty;% doubly-linked list of active verts
        RightNeighbor=vertex.empty;% doubly-linked list of active verts
        LeftF=face.empty;% face boardered by LeftV and ForwardV
        RightF=face.empty;% face boardered by RightV and ForwardV
        CollapsedF=face.empty;% face boardered by LeftV and RightV
        isbase=false;% Boolean - is this a base vertex (part of polygon)?
        status=1;% 1 - valid, 2 - behind left parent, 3 - behind right parent, 4 - invalid
        multi% pointer to itself if single, or to the associated multiVert
        ind=0;% index, used to construct vector output
    end
      
    methods 
        function xy=get.xy(thisV)
            % gets the multiVert position if there is one, otherwise
            % returns local xy, since multi points to itself for singles.
            xy=thisV.multi.xy;
        end
        
        function z=get.z(thisV)
            % gets the multiVert height if there is one, otherwise
            % returns local z, since multi points to itself for singles.
            z=thisV.multi.z;
        end
        
        function delta=get.delta(thisV)
            % gets the multiVert error if there is one, otherwise
            % returns local delta, since multi points to itself for singles.
            delta=thisV.multi.delta;
        end
        
        function thisV = vertex(LeftV,RightV,xy) 
            % constructor
            global scale
            thisV.multi=thisV;
            if nargin==3% base vertex, no parents
                thisV.isbase=true;
                thisV.xy=xy;
                thisV.z=0;
                thisV.delta=0;
                return
            end
            % find the intersection between two input vertices (xy is
            % ignored)
            thisV.delta=inf;
            thisV.LeftV=LeftV;
            thisV.RightV=RightV;
            
            thisV.LeftF=LeftV.LeftF;
            thisV.RightF=RightV.RightF;
            thisV.CollapsedF=LeftV.RightF;
            
            % don't intersect reflex base vertices
            if LeftV.isbase && RightV.isbase && ...
                    (LeftV.bi(1)*RightV.bi(2)-LeftV.bi(2)*RightV.bi(1))<0
                invalidate(thisV)
                return
            end
            biL=(thisV.LeftF.N+thisV.CollapsedF.N)*[0 -1;1 0];
            biR=(thisV.RightF.N+thisV.CollapsedF.N)*[0 -1;1 0];
            A=[thisV.LeftF.N 1
                thisV.RightF.N 1
                thisV.CollapsedF.N 1
                biL 0
                biR 0];
            b=[thisV.LeftF.N*thisV.LeftF.V.xy'
                thisV.RightF.N*thisV.RightF.V.xy'
                thisV.CollapsedF.N*thisV.CollapsedF.V.xy'
                biL*LeftV.xy'
                biR*RightV.xy'];
            condA=cond(A);
            if condA*eps>.1% too close to singular
                % invalid if inconsistent or pointed toward each other
                if rank(A)~=rank([A b]) || thisV.LeftV.bi*thisV.RightV.bi'<0
                    invalidate(thisV)
                    return
                else
                    % parallel - large error ensure multiVert is formed
                    thisV.delta=1/eps;% large, not inf
                end
            else
                % assume maximum error is 10 standard deviations
                thisV.delta=10*scale/6*eps*condA;
            end
            X=A\b;
            thisV.delta=thisV.delta+norm(A*X-b);% add residual error
            thisV.xy=X(1:2)';
            thisV.z=X(3);
            
            thisV.bi=bisectV(thisV);
            
            thisV.dL=(thisV.xy-LeftV.xy)*LeftV.bi';
            thisV.dR=(thisV.xy-RightV.xy)*RightV.bi';

            % the order only makes a difference when delta is large
            if (LeftV.xy-RightV.xy)*RightV.bi'<0
                BEHINDR=proximity(thisV,'dR',RightV);
                BEHINDL=proximity(thisV,'dL',LeftV);
            else
                BEHINDL=proximity(thisV,'dL',LeftV);
                BEHINDR=proximity(thisV,'dR',RightV);
            end
            if BEHINDL
                thisV.status=2;% behind left vertex
                thisV.z=inf;
            elseif BEHINDR
                thisV.status=3;% behind right vertex
                thisV.z=inf;
            end
        end
        
        function bisectBaseV(thisV)
            % Calculates the bisector vector for base vertices, which is
            % simple because the two faces adjoin, limiting the possible
            % angular arrangements.
            thisV.bi=-thisV.LeftF.N-thisV.RightF.N;
            thisV.bi=thisV.bi/norm(thisV.bi);
        end
        
        function bi=bisectV(thisV)
            % This is the most general and accurate calculation of a
            % bisector vector for a vertex, but requires knowledge of the
            % parent's bisectors, and hence is not suitable for base verts.
            bi=(thisV.LeftV.bi*thisV.RightV.bi')*...
                (-thisV.LeftF.N-thisV.RightF.N)+...
                (thisV.LeftV.bi(1)*thisV.RightV.bi(2)...
                -thisV.LeftV.bi(2)*thisV.RightV.bi(1))*...
                (thisV.RightF.N*[0 1;-1 0]-thisV.LeftF.N*[0 1;-1 0]);
            bi=bi/norm(bi);
        end
        
        function behind=proximity(thisV,d,prevV)
            % Determine if this intersection is behind its parents, and
            % therefore invalid. If it is close enough to a parent, it is
            % merged into it to form or add to a multiVert.
            behind=false;
            if thisV.(d)^2<=thisV.delta^2+prevV.delta^2% close enough to coincide
                if thisV.multi==thisV && prevV.multi==prevV% neither are multi-verts
                    thisMV=multiVert(thisV,prevV);
                    thisV.multi=thisMV;
                    prevV.multi=thisMV;
                else
                    if thisV.multi==thisV% thisV is not a multi-vert
                        thisMV=mergeMult(prevV.multi,thisV);
                    else
                        thisMV=mergeMult(thisV.multi,prevV);
                    end
                    [thisMV.Vlist.multi]=deal(thisMV);
                end
                thisV.(d)=0;
                if strcmp(d,'dL')
                    thisV.dR=(thisV.xy-thisV.RightV.xy)*thisV.RightV.bi';
                else
                    thisV.dL=(thisV.xy-thisV.LeftV.xy)*thisV.LeftV.bi';
                end
            elseif thisV.(d)<0
                behind=true;
            end
        end
        
        function linkV(LeftV,RightV)
            LeftV.RightV=RightV;
            RightV.LeftV=LeftV;
            LeftV.RightNeighbor=RightV;
            RightV.LeftNeighbor=LeftV;    
        end
        
        function linkVF(LeftV,thisF)
            LeftV.RightF=thisF;
            LeftV.RightV.LeftF=thisF;
        end
               
        function finalLink(V1,V2)
            V1.ForwardV=V2;
            V2.ForwardV=V1;
        end
        
        function addV(thisV)
            thisV.LeftV.ForwardV=thisV;
            thisV.RightV.ForwardV=thisV;
            thisV.LeftV.LeftNeighbor.RightNeighbor=thisV;
            thisV.RightV.RightNeighbor.LeftNeighbor=thisV;
            thisV.LeftNeighbor=thisV.LeftV.LeftNeighbor;
            thisV.RightNeighbor=thisV.RightV.RightNeighbor;
        end
        
        function removeV(badV)
            left=badV.LeftV;
            right=badV.RightV;
            left.RightNeighbor=right;
            right.LeftNeighbor=left;
            left.LeftNeighbor=badV.LeftNeighbor;
            right.RightNeighbor=badV.RightNeighbor;
            badV.LeftNeighbor.RightNeighbor=left;
            badV.RightNeighbor.LeftNeighbor=right;
            badV.status=4;
        end
        
        function flagRemoval(thisV,badV)
            thisV.toRemove=badV;
        end
        
        function invalidate(thisV)
            thisV.dL=inf;
            thisV.dR=inf;
            thisV.status=4;
        end       
        
        function V=keyhole(lowerV,bottomV,bottomP)
            % create new vertices
            V(1)=vertex([],[],bottomP);
            V(2)=vertex([],[],bottomV.xy);
            V(3)=vertex([],[],bottomP);
            % create new faces
            F(1)=face(V(1),[1 0]);
            F(2)=face(V(2),[-1 0]);   
            % link vertices
            linkV(bottomV.LeftV,V(2));
            linkV(V(3),lowerV.RightV);
            linkV(lowerV,V(1));
            linkV(V(1),bottomV);            
            linkV(V(2),V(3));
            % link faces
            linkVF(V(1),F(1));
            linkVF(V(2),F(2));
            V(1).LeftF=V(1).LeftV.RightF;
            V(2).LeftF=V(2).LeftV.RightF;
            V(3).RightF=V(3).RightV.LeftF;
            % calculate bisectors
            bisectBaseV(V(1));
            bisectBaseV(V(2));
            bisectBaseV(V(3));            
            bisectBaseV(bottomV);
        end
        
        function unKeyhole(bottomV,V)
            % remove V(1)         
            V(3).LeftF=V(1).LeftF;            
            linkV(V(1).LeftV,V(3)); 
            % remove V(2)
            bottomV.LeftF=V(2).LeftF;                
            linkV(V(2).LeftV,bottomV);
            % recalculate bisector angle           
            bisectBaseV(bottomV);
            % open keyhole
            removeFaceV(V(1),bottomV);
            removeFaceV(V(2),V(3));
            % connect top
            first=nextV(V(2).ForwardV,V(2),-1);
            bottomV.LeftNeighbor=first;
            first.RightNeighbor=bottomV;
            bottomV.ForwardV=[];
            % connect bottom and remove V(3)
            linkV(V(3).LeftV,V(3).RightV);
            first=nextV(V(1).ForwardV,V(1),-1);
            first.RightNeighbor=V(3).RightNeighbor;
            V(3).RightNeighbor.LeftNeighbor=first;
            invalidate(V(3));
        end
        
        function removeFaceV(startV,endV)
            % Remove all vertices around a face between startV and endV
            last=startV;
            current=startV.ForwardV;
            currentOut=nextV(current,last,-1);
            permuteV(currentOut,current);
            currentOut.ForwardV=[];
            invalidate(startV);
            while current~=endV;   
                next=nextV(current,last,1);  
                currentOut=nextV(current,last,-1);
                invalidate(current);               
                if next~=endV
                    nextOut=nextV(next,current,-1);
                    permuteV(nextOut,next);
                    nextOut.ForwardV=[];
                    currentOut.LeftNeighbor=nextOut;
                    nextOut.RightNeighbor=currentOut;
                else
                    currentOut.LeftNeighbor=next;
                    next.RightNeighbor=currentOut;
                end
                last=current;
                current=next;
            end
        end
        
        function permuteV(thisV,forwardV)
            % Rotate the vertex so that the ForwardV property points
            % toward the input vertex, forwardV.
            thisV.QL=vertex.empty;
            thisV.QR=vertex.empty;
            LFR=[thisV.LeftV thisV.ForwardV thisV.RightV];
            i=find(LFR==forwardV);
            if i~=2% not already correct
                j=[i:3 1:i-1];
                thisV.ForwardV=LFR(j(1));
                thisV.RightV=LFR(j(2));
                thisV.LeftV=LFR(j(3));
                
                RCLF=[thisV.RightF;thisV.CollapsedF;thisV.LeftF];
                thisV.CollapsedF=RCLF(j(1));
                thisV.LeftF=RCLF(j(2));
                thisV.RightF=RCLF(j(3));
                
                thisV.bi=thisV.RightF.N*[0 1;-1 0]-thisV.LeftF.N*[0 1;-1 0];
                thisV.bi=thisV.bi/norm(thisV.bi);
                thisV.dL=norm(thisV.xy-thisV.LeftV.xy);
                thisV.dR=norm(thisV.xy-thisV.RightV.xy);
            end
        end
        
        function V=nextV(thisV,lastV,k)
            % RightV, LeftV, ForwardV are always in clockwise order. nextV
            % uses this to find the next vertex k positions from lastV
            % around thisV. Positive k is clockwise. 
            V=[thisV.LeftV thisV.RightV thisV.ForwardV];
            i=find(V==lastV,1);
            V=V(mod(i-1+k,3)+1);
        end
        
        function addQ(thisQ)
            thisQ.LeftV.QL(end+1)=thisQ;
            thisQ.RightV.QR(end+1)=thisQ;
        end
        
        function removeQ(thisQ)
            thisQ.LeftV.QL(thisQ.LeftV.QL==thisQ)=[];
            thisQ.RightV.QR(thisQ.RightV.QR==thisQ)=[];
        end
        
        function checkQ(thisV)
            % removes candidate vertices that are not valid due to not
            % being connected in the active neighbor list
            keep=true(size(thisV.QL));
            for i=1:length(thisV.QL)
                keep(i)=isConnected(thisV.QL(i));
            end
            thisV.QL=thisV.QL(keep);
            keep=true(size(thisV.QR));
            for i=1:length(thisV.QR)
                keep(i)=isConnected(thisV.QR(i));
            end
            thisV.QR=thisV.QR(keep);
        end
        
        function bool=isConnected(thisQ)
            % returns true if the input vertex is connected in the active
            % neighbor list
            LV=thisQ.LeftV;
            RV=thisQ.RightV;
            if ~isempty(thisQ.toRemove)% use first removal instead of parent
                if thisQ.toRemove(1).RightF.V==thisQ.CollapsedF.V
                    LV=thisQ.toRemove(1);
                else
                    RV=thisQ.toRemove(1);
                end
            end
            bool=LV.RightNeighbor==RV && RV.LeftNeighbor==LV;
        end
        
        function bool=isFirst(thisQ)
            % returns true if thisQ is the nearest candidate of its parents
            checkQ(thisQ.LeftV);% remove stale candidates
            checkQ(thisQ.RightV);
            Lmin=min([thisQ.LeftV.QL.dL thisQ.LeftV.QR.dR]+...
                [thisQ.LeftV.QL.delta thisQ.LeftV.QR.delta]);
            Rmin=min([thisQ.RightV.QL.dL thisQ.RightV.QR.dR]+...
                [thisQ.RightV.QL.delta thisQ.RightV.QR.delta]);
            bool=(isempty(Lmin) || thisQ.dL-thisQ.delta<Lmin) && ...
                (isempty(Rmin) || thisQ.dR-thisQ.delta<Rmin);
        end
        
        function index(V,i)
            V.ind=i;
        end
        
        function showNeighbors(startV)
            % Function to plot the current loop of active vertices,
            % following the doubly-linked list of neighbors.
            list=[startV startV.RightNeighbor];
            while list(end)~=startV
                list(end+1)=list(end).RightNeighbor;
            end
            xy=vertcat(list.xy);
            plot(xy(:,1),xy(:,2),'rx:')
        end
        
        function plot(V,varargin)
            % Plot partial skeleton to figure 1. Takes the same extra 
            % arguments as a normal plot command.
            global bb
            figure(1)
            clf
            hold on
            X=zeros(0,2);
            Y=zeros(0,2);
            for i=1:length(V)
                if V(i).status==1
                    X=[X
                        V(i).xy(1) V(i).LeftV.xy(1)
                        V(i).xy(1) V(i).RightV.xy(1)];
                    Y=[Y
                        V(i).xy(2) V(i).LeftV.xy(2)
                        V(i).xy(2) V(i).RightV.xy(2)];
                    if ~isempty(V(i).ForwardV)
                        X=[X
                            V(i).xy(1) V(i).ForwardV.xy(1)];
                        Y=[Y
                            V(i).xy(2) V(i).ForwardV.xy(2)];
                    end
                end
            end           
            plot(X',Y',varargin{:})
            axis equal
            axis(bb)
%             axis off
        end
    end
end
