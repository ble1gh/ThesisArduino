clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 340;
theta = 0;

%read in data
Rc_tab = readtable('Rc_range.csv');
Rc_data = -table2array(Rc_tab);
predicted_tab = readtable('predicted_theta_0.00.csv');
predicted_data = table2array(predicted_tab);
output_tab = readtable('output_theta_0.00.csv');
output_data = table2array(output_tab);
pos_data = output_data(:,7:9)';

%calculate error
model_error = norm(pos_data-predicted_data)

%curve fit
[xfit_Rc, gofx] = fit(pos_data(1,:)',Rc_data','exp2');

figure(1)
plot(xfit_Rc,pos_data(1,:),Rc_data,'.')
hold on; grid on;
plot(predicted_data(1,:),Rc_data)
legend('Data','Curve Fit','Prediction','Location','east')
xlabel('x distance from center (mm)')
ylabel('Rc')
title('Rc vs x, Prediction vs Reality')

%generate new points to test
R_fit = [xfit_Rc.a; xfit_Rc.b; xfit_Rc.c; xfit_Rc.d];
x = linspace(3,20,20);
R = R_fit(1)*exp(R_fit(2)*x) + R_fit(3)*exp(R_fit(4)*x);

%initialize prediction and input vars
pred = zeros(3,length(R));
dl = zeros(4,length(R)); 
u = ones(4,length(dl))*setmid;

%predicted xz
pred(1,:) = -cos(0)*R.*(1-cos(L./R));
pred(2,:) = -sin(0)*R.*(1-cos(L./R));
pred(3,:) = R.*sin(L./R);

%cable length change
dl(1,:) = (L./R.*(R-d*cos(theta))-L);
dl(2,:) = (L./R.*(R+d*cos(theta))-L);
dl(3,:) = (L./R.*(R-d*sin(theta))-L);
dl(4,:) = (L./R.*(R+d*sin(theta))-L);

%convert to motor values
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u = u-dPWM;

results_tab = readtable("transformed_output_0516_xz.csv");
results_data = table2array(results_tab);

figure(2)
plot(results_data(1,3:end),R(3:end),'.',pred(1,3:end),R(3:end),predicted_data(1,:),Rc_data)

figure(1)
plot(results_data(1,3:end),R(3:end),'g.','MarkerSize',20)

figure(3)
plot3(results_data(1,:),results_data(2,:),results_data(3,:),'.',pred(1,:),pred(2,:),pred(3,:),'.')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
grid on; axis equal;

results_error = norm(results_data(1,:)-pred(1,:))