function []=plotslice(s,z)

hold on
color=get(0,'defaultaxescolororder');
for i=1:length(s)
    for j=1:length(s(i).loop)
        xyz=[s(i).loop(j).xy zeros(size(s(i).loop(j).xy,1),1)+z(i)
            s(i).loop(j).xy(1,:) z(i)];% close loop
        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'color',color(mod(j-1,7)+1,:));
        plot3(xyz(1,1),xyz(1,2),xyz(1,3),'.');
        plot3(xyz(2,1),xyz(2,2),xyz(2,3),'g.');
%         fill3(xyz(:,1),xyz(:,2),xyz(:,3),'r','facealpha',0.2);
    end
end
axis equal
axis vis3d
view(3);