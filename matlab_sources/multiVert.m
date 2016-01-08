classdef multiVert < handle
    
   properties(SetAccess = private)
       xy% horizontal position
       z% distance from associated edges (height of roof)
       delta% uncertainty, reduced by averaging vertices
       Vlist% list of vertices pointing to this multiVert
   end
   
   methods
       function thisMV = multiVert(V1,V2)
           setVals(thisMV,V1,V2);
           thisMV.Vlist=[V1 V2];
       end
       
       function thisMV=mergeMult(thisMV,V)         
           if V.multi==V% if V is not a multiVert
               Vs=V;
           else
               Vs=V.multi.Vlist;
           end
           thisMV.Vlist=unique([thisMV.Vlist Vs]);
           setVals(thisMV,thisMV,V);
       end
       
       function setVals(thisMV,V1,V2)
           a=V2.delta^2/(V1.delta^2+V2.delta^2);% convex parameter
           a=min(max(a,0),1);% clamp to bounds
           thisMV.xy=a*V1.xy+(1-a)*V2.xy;% error-weighted average
           thisMV.z=a*V1.z+(1-a)*V2.z;
           thisMV.delta=a*V1.delta+(1-a)*V2.delta;% improved uncertainty
       end
   end
end