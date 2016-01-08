function poly=rotatePoly(poly,R)
for i=1:length(poly)
    poly(i).xy=poly(i).xy*R;
end