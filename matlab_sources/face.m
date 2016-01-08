classdef face < handle
    
    properties(SetAccess = private)
        V=vertex.empty;% Left base vertex of this face
        N% 2D unit normal vector pointing away from the polygon
    end
      
    methods 
        function thisF=face(LV,normal)
            thisF.V=LV;
            thisF.N=normal;
        end
    end
end