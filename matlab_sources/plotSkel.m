function plotSkel(skel,varargin)

if ~isempty(skel.xy)
    halfedge=checkSkel(skel);
    x=reshape(skel.xy(halfedge,1),[],2);
    y=reshape(skel.xy(halfedge,2),[],2);
    
    plot(x',y',varargin{:})
    axis equal
    axis off
    hold on
end