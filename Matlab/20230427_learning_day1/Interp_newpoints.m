close all; clear;

%load data into var
history_tab = readtable('output_input_history_0427.csv','Delimiter',',');
ustar = table2array(history_tab);

% here's a query point (used an existing point to test)
xq = -20;
yq = -5;

% find associated z value for (x,y) pair
x = sortrows(ustar(1,:)')';
y = sortrows(ustar(2,:)')';
[X,Y] = meshgrid(x,y);
Z = griddata(ustar(1,:),ustar(2,:),ustar(3,:),X,Y);
zq = interp2(X,Y,Z,xq,yq,"linear")

%create interpolant objects for each motor given (x,y,z) position
F1 = scatteredInterpolant(ustar(1:3,:)',ustar(4,:)','linear','linear');
F2 = scatteredInterpolant(ustar(1:3,:)',ustar(5,:)','linear','linear');
F3 = scatteredInterpolant(ustar(1:3,:)',ustar(6,:)','linear','linear');
F4 = scatteredInterpolant(ustar(1:3,:)',ustar(7,:)','linear','linear');

%evaluate the interpolant objects to find motor values given a (x,y,z) pos
% m1 = round(F1(xq,yq,zq))
% m2 = round(F2(xq,yq,zq))
% m3 = round(F3(xq,yq,zq))
% m4 = round(F4(xq,yq,zq))

m1 = F1(xq,yq,zq)
m2 = F2(xq,yq,zq)
m3 = F3(xq,yq,zq)
m4 = F4(xq,yq,zq)


%interp3(X,Y,Z,)

figure(1)
plot3(X,Y,Z,'.','Color',"#0072BD")
hold on; plot3(xq,yq,zq,'.g','MarkerSize',20); hold off;
title('3D Position Data Used for Interpolation')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;