%evaluation of accuracy of curve-fit control method
clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;

%set frequency of sampling
res_curve = 18;
res_theta = 40;
iterations = 3;

points_record_tab = readtable("points_record_curvefit_test.csv");
points_record = table2array(points_record_tab);
requests_tab = readtable("points_request_curvefit_test.csv");
requests = table2array(requests_tab);
zfits_tab = readtable('zfits_poly3_full_circle.csv');
zfits = table2array(zfits_tab);

theta = linspace(0,2*pi,res_theta);
r_req = linspace(0,30,res_curve);

fit_theta = zfits(1,:);
az = interpn(fit_theta,zfits(2,:),theta);
bz = interpn(fit_theta,zfits(3,:),theta);
cz = interpn(fit_theta,zfits(4,:),theta);
dz = interpn(fit_theta,zfits(5,:),theta);

pos_req = zeros(3,res_curve*res_theta);
for i = 1:res_theta
    pos_req(:,(i-1)*res_curve+1:i*res_curve) = [r_req.*cos(theta(i)); r_req.*sin(theta(i)); 
    az(i)*r_req.^3 + bz(i)*r_req.^2 + cz(i)*r_req.^1 + dz(i)];
end

%calculate error
error = vecnorm(points_record(1:3,:)-pos_req);
avg_error = mean(error)
[max_error, max_idx] = max(error)
[min_error, min_idx] = min(error)

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

figure(2)
%plot3(requests(1,:),requests(2,:),requests(3,:),'.')
plot3(points_record(1,:),points_record(2,:),points_record(3,:),'.')
hold on; grid on;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Data-based control output data','FontSize',18)
axis equal

figure(3)
plot(requests(1,1:18),requests(3,1:18),'.')
hold on; grid on;
plot(points_record(1,end-18:end),points_record(3,end-18:end),'.')

%sort data into gridded form
x = sortrows(pos_req(1,:)')';
y = sortrows(pos_req(2,:)')';
[X,Y] = meshgrid(x,y);
Z = griddata(pos_req(1,:),pos_req(2,:),pos_req(3,:),X,Y);
E = griddata(pos_req(1,:),pos_req(2,:),error,X,Y);

thetamistake = 2*pi+atan2(pos_req(2,595:600),pos_req(1,595:600))

figure(4)
surf(X,Y,Z,E,'EdgeColor','none')
colorbar
hold on
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
plot3(pos_req(1,:),pos_req(2,:),pos_req(3,:),'.')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Data-based control surface with error (mm) as color','FontSize',18)
axis equal
clim([0 14]);

figure(5)
plot(points_record(1,:),points_record(2,:),'.')
hold on; grid on;
plot(requests(1,:),requests(2,:),'.')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Data-based control output data','FontSize',18)
axis equal
legend('data','prediction')