%Data analysis from today's data run, loads data from files, calculates
%error between requested and resulting, and then meshes a surface and plots
%area reached and corresponding error

close all; clear;

% load transformed xyz data and requested xyz points
points_tab = readtable('transformed_output_0512.csv','Delimiter',',');
spinedata_in_modelspace = table2array(points_tab);
predicted_tab = readtable('predicted.csv','Delimiter',',');
requested_points = table2array(predicted_tab);
neutralpoint = [0 0 56]';
requested_points = [requested_points neutralpoint];
motor_tab = readtable('input.csv','Delimiter',',');
u = table2array(motor_tab);
neutral_motor = [300 300 300 300]';
u = [u neutral_motor];
output_input_history = [spinedata_in_modelspace; u];
writematrix(output_input_history,'output_input_history_0512.csv')

%error = sqrt((spinedata_in_modelspace(1,:)-requested_points(1,:)).^2+(spinedata_in_modelspace(2,:)-requested_points(2,:)).^2+(spinedata_in_modelspace(3,:)-requested_points(3,:)).^2);
error = vecnorm(spinedata_in_modelspace(1:3,:)-requested_points(1:3,:));

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);  % not sure why we can pass an array of more than 3 cols for vertices here, but it is OK with homogeneous coords...


%convert output data to grid for surf
[Xsurf,Ysurf] = meshgrid(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:));
Zsurf = griddata(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),Xsurf,Ysurf);
Esurf = griddata(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),error,Xsurf,Ysurf);

figure(1)
%surf(Xsurf,Ysurf,Zsurf,Esurf,'FaceColor',"#4DBEEE",'FaceAlpha',0.5,'EdgeColor','none')
surf(Xsurf,Ysurf,Zsurf,Esurf,'EdgeColor','none')
colorbar
hold on
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','MarkerSize',10)
%plot3(requested_points(1,:),requested_points(2,:),requested_points(3,:),'.','MarkerSize',10)
title('Model-Based Control with Error in mm from Expected Position')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;
hold off

%for debugging
point = 2;
figure(2)
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','MarkerSize',10)
hold on;
plot3(requested_points(1,:),requested_points(2,:),requested_points(3,:),'.','MarkerSize',10)
%plot3(spinedata_in_modelspace(1,point),spinedata_in_modelspace(2,point),spinedata_in_modelspace(3,point),'.m',requested_points(1,point),requested_points(2,point),requested_points(3,point),'.g','MarkerSize',20)
legend('data','model')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;
hold off