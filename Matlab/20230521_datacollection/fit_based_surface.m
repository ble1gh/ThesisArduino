clear; close all;

x_coeffs_tab = readtable("x_coeffs_full.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs_tab = readtable("zfits_poly3_full_circle.csv");
z_coeffs = table2array(z_coeffs_tab);

res_curve = 18;
r_req = linspace(0,30,res_curve);
theta = x_coeffs(1,:);
res_theta = length(theta);
pos = zeros(3,res_theta*res_curve);

fit_theta = z_coeffs(1,:);
az = z_coeffs(2,:);
bz = z_coeffs(3,:);
cz = z_coeffs(4,:);
dz = z_coeffs(5,:);

for i = 1:res_theta
    pos(:,(i-1)*res_curve+1:i*res_curve) = [r_req*cos(theta(i)); r_req*sin(theta(i));
        az(i)*r_req.^3 + bz(i)*r_req.^2 + cz(i)*r_req.^1 + dz(i)];
end

%sort data into arranged (X,Y) grid with associated Z values
[z_unique, p_unique] = groupsummary(pos(3,:)',pos(1:2,:)',"mean");
xy = [p_unique{:}];
x = xy(:,1);
y = xy(:,2);
[X,Y] = meshgrid(x,y);
Z = griddata(x,y,z_unique,X,Y);

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

figure(1)
surf(X,Y,Z)
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);

