close all; clear;

%load data into var
history_tab = readtable('points_record_circle.csv','Delimiter',',');
ustar = table2array(history_tab);

% here's a query point (used an existing point to test)
xq = -20;
yq = -5;

% find associated z value for (x,y) pair
[z_unique, p_unique] = groupsummary(ustar(3,:)',ustar(1:2,:)',"mean");
xy = [p_unique{:}];
x = xy(:,1);
y = xy(:,2);
xsorted = sortrows(x)';
ysorted = sortrows(y)';
[Xdata,Ydata] = meshgrid(xsorted,ysorted);
Zdata = griddata(x,y,z_unique,Xdata,Ydata);

zq1 = interp2(Xdata,Ydata,Zdata,xq,yq,"linear")

x_coeffs_tab = readtable("x_coeffs_full.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs_tab = readtable("zfits_averaged.csv");
z_coeffs = table2array(z_coeffs_tab);

res_curve = 18;
theta = x_coeffs(1,:);
r_req = linspace(0,30,res_curve);
res_theta = length(theta);
pos = zeros(3,res_theta*res_curve);

for i = 1:res_theta
    pos(:,(i-1)*res_curve+1:i*res_curve) = [r_req*cos(theta(i)); r_req*sin(theta(i));
        z_coeffs(2,i)*r_req.^4 + z_coeffs(3,i)*r_req.^3 + z_coeffs(4,i)*r_req.^2 + z_coeffs(5,i)*r_req.^1 + z_coeffs(6,i)];
end

%sort data into arranged (X,Y) grid with associated Z values
[z_unique, p_unique] = groupsummary(pos(3,:)',pos(1:2,:)',"mean");
xy = [p_unique{:}];
x = xy(:,1);
y = xy(:,2);
zcurve_data = interp2(Xdata,Ydata,Zdata,x,y,"linear");
xsorted = unique(x)';
ysorted = unique(y)';
[Xcurve,Ycurve] = meshgrid(xsorted,ysorted);
Zcurve = griddata(x,y,z_unique,Xcurve,Ycurve);
Zcurve_data = griddata(x,y,zcurve_data,Xcurve,Ycurve);

curve_fit_z_error = abs(Zcurve_data-Zcurve);

zq2 = interp2(Xcurve,Ycurve,Zcurve,xq,yq,"natural")

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

%%

% m1 = F1(xq,yq,zq);
% m2 = F2(xq,yq,zq);
% m3 = F3(xq,yq,zq);
% m4 = F4(xq,yq,zq);


% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

figure(1)
%surf(Xcurve,Ycurve,Zcurve,'EdgeColor','none','FaceColor',"#4DBEEE",'FaceAlpha',0.7)
surf(Xcurve,Ycurve,Zcurve,curve_fit_z_error,'EdgeColor','none')
colorbar
hold on; 
plot3(xq,yq,zq2,'.g','MarkerSize',20);
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
title('3D Position Data Used for Interpolation')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;


