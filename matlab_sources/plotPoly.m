function plotPoly(poly,varargin)
hold on
for i=1:length(poly)
    plot(poly(i).xy([1:end,1],1),poly(i).xy([1:end,1],2),varargin{:})
end