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
zfits_tab = readtable('zfits_3runs.csv');
zfits = table2array(zfits_tab);
az = zfits(2,:);
bz = zfits(3,:);
cz = zfits(4,:);
dz = zfits(5,:);
ez = zfits(6,:);

theta = linspace(0,2*pi,res_theta);
r_req = linspace(0,30,res_curve);
pos_req = zeros(3,res_curve*res_theta);
for i = 1:res_theta
    pos_req(:,(i-1)*res_curve+1:i*res_curve) = [r_req.*cos(theta(i)); r_req.*sin(theta(i)); 
    az(i)*r_req.^4 + bz(i)*r_req.^3 + cz(i)*r_req.^2 + dz(i)*r_req.^1 + ez(i)];
end

%calculate error
error = vecnorm(points_record(1:3,:)-pos_req);
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
plot3(pos_req(1,:),pos_req(2,:),pos_req(3,:),'.')
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);

figure(2)
plot3(pos_req(1,:),pos_req(2,:),pos_req(3,:),'.')
hold on; grid on;
plot3(points_record(1,:),points_record(2,:),points_record(3,:),'.')
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
legend('Requests','Position Data')

figure(3)
plot3(pos_req(1,:), pos_req(2,:), pos_req(3,:),'.')