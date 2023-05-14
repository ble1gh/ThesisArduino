close all; clear;

% load transformed xyz data
points_tab = readtable('transformed_output_0512.csv','Delimiter',',');
spinedata_in_modelspace = table2array(points_tab);
predicted_tab = readtable('predicted.csv','Delimiter',',');
requested_points = table2array(predicted_tab);

[Xsurf,Ysurf] = meshgrid(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:));
Zsurf = griddata(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),Xsurf,Ysurf);

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);  % not sure why we can pass an array of more than 3 cols for vertices here, but it is OK with homogeneous coords...


%plot
figure(1)
%nicer green:
%plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','Color',"#77AC30");
%easier to see:
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'r.','MarkerSize',20);
hold on; 
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
surf(Xsurf,Ysurf,Zsurf,'FaceColor',"#4DBEEE",'FaceAlpha',0.5,'EdgeColor','none')
axis equal; hold on; grid on;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x')
ylabel('y')
zlabel('z')
hold off;

%simulated reachable area
L = 56; % (m) from spine architecture 

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(3*L/(2*pi));
maxR = log10(200);
Rc = [logspace(minR,maxR,50)];


theta = linspace(0,2*pi,50);
phi = L./Rc;

%initiate position variables
x = zeros(length(theta),length(Rc));
y = zeros(length(theta),length(Rc));
z = zeros(length(theta),length(Rc));

%for a given theta, an arc is defined along the range of possible bend
%radiuses which each give an (x,y,z) coordinate
for i = 1:length(theta)
    x(i,:) = cos(theta(i))*Rc.*(1-cos(L./Rc));
    y(i,:) = sin(theta(i))*Rc.*(1-cos(L./Rc));
    z(i,:) = Rc.*sin(L./Rc);
end
% FV.vertices = x';
% points = spinedata_in_modelspace';
% [distances,surface_points] = point2trimesh(FV, 'QueryPoints', points);

%expectedZ = interpn(x,y,z,spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),'linear',0);


%plot result
figure(1)
%surf(x,y,z,'FaceColor','interp','FaceAlpha',0.9,'EdgeColor','none')
surf(x,y,z,'FaceColor',"#4DBEEE",'FaceAlpha',0.7)
%plot3(X,Y,Z,'.','Color',"#4DBEEE")
title('Simulated vs Actual Reachable Region')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;

hold off;

% figure(3)
% surf(x,y,z,'FaceColor',"#4DBEEE",'FaceAlpha',0.7)
% title('Simulated Reachable Region, 3D')
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')
