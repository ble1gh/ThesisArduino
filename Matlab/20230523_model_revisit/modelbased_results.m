clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;

points_record_tab = readtable("points_record_model_test.csv");
points_record = table2array(points_record_tab);
requests_tab = readtable("points_request_model_test.csv");
requests = table2array(requests_tab);

%calculate error
error = vecnorm(points_record(1:3,:)-requests);
avg_error = mean(error)

%sort data into gridded form
x = sortrows(points_record(1,:)')';
y = sortrows(points_record(2,:)')';
[X,Y] = meshgrid(x,y);
Z = griddata(points_record(1,:),points_record(2,:),points_record(3,:),X,Y);
E = griddata(points_record(1,:),points_record(2,:),error,X,Y);

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

figure(1)
surf(X,Y,Z,E,'EdgeColor','none')
colorbar
hold on
plot3(requests(1,:),requests(2,:),requests(3,:),'.')
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);

%sort data into gridded form
x = sortrows(requests(1,:)')';
y = sortrows(requests(2,:)')';
[X,Y] = meshgrid(x,y);
Z = griddata(requests(1,:),requests(2,:),requests(3,:),X,Y);
E = griddata(requests(1,:),requests(2,:),error,X,Y);

figure(2)
surf(X,Y,Z,E,'EdgeColor','none')
colorbar
hold on; grid on;
plot3(requests(1,:),requests(2,:),requests(3,:),'.')
%plot3(points_record(1,:),points_record(2,:),points_record(3,:),'.')
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('CC-based control surface with error (mm) as color','FontSize',18)
axis equal
clim([0 14])

figure(3)
plot3(points_record(1,:),points_record(2,:),points_record(3,:),'.')
hold on; grid on;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('CC-based control data','FontSize',18)
axis equal