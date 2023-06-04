clear; close all;

x_coeffs_tab = readtable("x_coeffs_full.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs_tab = readtable("zfits_poly3_full_circle.csv");
z_coeffs = table2array(z_coeffs_tab);

res_curve = 50;
r_req = linspace(0,30,res_curve);
res_theta = 150;
theta = linspace(0,2*pi,res_theta);
pos = zeros(3,res_theta*res_curve);

fit_theta = z_coeffs(1,:);
az = interpn(fit_theta,z_coeffs(2,:),theta);
bz = interpn(fit_theta,z_coeffs(3,:),theta);
cz = interpn(fit_theta,z_coeffs(4,:),theta);
dz = interpn(fit_theta,z_coeffs(5,:),theta);

for i = 1:res_theta
    pos(:,(i-1)*res_curve+1:i*res_curve) = [r_req*cos(theta(i)); r_req*sin(theta(i));
        az(i)*r_req.^3 + bz(i)*r_req.^2 + cz(i)*r_req.^1 + dz(i)];
end

%sort data into arranged (X,Y) grid with associated Z values
[z_unique, p_unique] = groupsummary(pos(3,:)',pos(1:2,:)',"mean");
xy = [p_unique{:}];
x = xy(:,1);
y = xy(:,2);
workspace = [x'; y'; z_unique'];

writematrix(workspace,"known_surface.csv")

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

figure(1)
plot3(workspace(1,:),workspace(2,:),workspace(3,:),'.')
grid on;

pen = [20, -20, 70]';

dist = vecnorm(workspace-pen);
[val,indx] = min(dist);
req = workspace(:,indx)


figure(1)
hold on
plot3(pen(1),pen(2),pen(3),'g.','MarkerSize',30);
plot3(req(1),req(2),req(3),'m.','MarkerSize',30)
axis equal
title('Discretized workspace with example pen-position and objective');
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
legend('reachable points','pen','closest point')