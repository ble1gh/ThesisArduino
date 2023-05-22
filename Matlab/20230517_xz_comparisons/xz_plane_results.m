clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;
theta = 0;
res = 18;

%read in data
Rc_tab = readtable('Rc_data.csv');
Rc_data = -table2array(Rc_tab);
predicted_tab = readtable('modelb_predicted.csv');
predicted_data = table2array(predicted_tab);
output_tab = readtable('transformed_output_0517_xz.csv');
output_data = table2array(output_tab);
motorval_tab = readtable('modelb_input.csv');
motorvals = table2array(motorval_tab);
dl_tab = readtable('modelb_cable_dl.csv');
dl = table2array(dl_tab);

%predicted data has no y coordinates
predicted_data = [predicted_data(1,:); zeros(1,length(predicted_data));predicted_data(2,:)];

figure(1)
plot(output_data(1,:),output_data(3,:),'.')
hold on; grid on; axis equal;
plot(predicted_data(1,:),predicted_data(3,:),'.')
xlabel('x (mm)')
ylabel('y (mm)')
title('Model-Based Sweep, Predicted vs Data')
legend('Data','Prediction')

figure(2)
plot(output_data(1,:),dl(1,:),'.')
hold on; grid on;
plot(output_data(1,:),dl(2,:),'.')
xlabel('x (mm)')
ylabel('cable length change (mm)')
title('Model-Based Sweep, x Displacement vs Cable Length Change')
legend('Pulling motor','Giving motor')


[xfit_motor1, gofx] = fit(output_data(1,:)',dl(1,:)','poly4')

figure(3)
plot(xfit_motor1,output_data(1,:),dl(1,:))
grid on;
xlabel('x (mm)')
ylabel('dl (mm)')
title('In Plane x Displacement vs Cable Length change')

k = 1;
fit_coeffs = [xfit_motor1.p1; xfit_motor1.p2; xfit_motor1.p3; xfit_motor1.p4; xfit_motor1.p5];

writematrix(fit_coeffs,'fitted_coeffs_v1.csv')


pred_fitv1_tab = readtable('fitted_v1_predicted.csv');
pred_fitv1_data = table2array(pred_fitv1_tab);
output_fitv1_tab = readtable('fitted_v1_output_transformed.csv');
output_fitv1_data = table2array(output_fitv1_tab);
motorval_fitv1_tab = readtable('fitted_v1_input.csv');
motorvals_fitv1 = table2array(motorval_fitv1_tab);
dl_fitv1_tab = readtable('fitted_v1_cable_dl.csv');
dl_fitv1 = table2array(dl_fitv1_tab);

figure(4)
plot(output_fitv1_data(1,:),output_fitv1_data(3,:))

[xfit_motor1, gofx] = fit(output_fitv1_data(1,:)',dl_fitv1(1,:)','poly4')

k = 2;
fit_coeffs(:,k) = [xfit_motor1.p1; xfit_motor1.p2; xfit_motor1.p3; xfit_motor1.p4; xfit_motor1.p5];

writematrix(fit_coeffs,'fitted_coeffs_v2.csv')

figure(3)
hold on
plot(xfit_motor1,output_fitv1_data(1,:),dl_fitv1,'.')

%iteration 3
pred_fitv2_tab = readtable('fitted_v2_predicted.csv');
pred_fitv2_data = table2array(pred_fitv2_tab);
output_fitv2_tab = readtable('fitted_v2_output_transformed.csv');
output_fitv2_data = table2array(output_fitv2_tab);
motorval_fitv2_tab = readtable('fitted_v2_input.csv');
motorvals_fitv2 = table2array(motorval_fitv2_tab);
dl_fitv2_tab = readtable('fitted_v2_cable_dl.csv');
dl_fitv2 = table2array(dl_fitv2_tab);

figure(4)
plot(output_fitv2_data(1,:),output_fitv2_data(3,:))

[xfit_motor1, gofx] = fit(output_fitv2_data(1,:)',dl_fitv2(1,:)','poly4')

k = 3;
fit_coeffs(:,k) = [xfit_motor1.p1; xfit_motor1.p2; xfit_motor1.p3; xfit_motor1.p4; xfit_motor1.p5];

writematrix(fit_coeffs,'fitted_coeffs_v3.csv')

figure(3)
hold on
plot(xfit_motor1,output_fitv2_data(1,:),dl_fitv2,'.')

%iteration 4
pred_fitv3_tab = readtable('fitted_v3_predicted.csv');
pred_fitv3_data = table2array(pred_fitv3_tab);
output_fitv3_tab = readtable('fitted_v3_output_transformed.csv');
output_fitv3_data = table2array(output_fitv3_tab);
motorval_fitv3_tab = readtable('fitted_v3_input.csv');
motorvals_fitv3 = table2array(motorval_fitv3_tab);
dl_fitv3_tab = readtable('fitted_v3_cable_dl.csv');
dl_fitv3 = table2array(dl_fitv3_tab);

figure(4)
plot(output_fitv3_data(1,:),output_fitv3_data(3,:))

[xfit_motor1, gofx] = fit(output_fitv3_data(1,:)',dl_fitv3(1,:)','poly4')
[zfit_motor1, gofz] = fit(output_fitv2_data(1,:)',output_fitv2_data(3,:)','poly4')

k = 4;
fit_coeffs(:,k) = [xfit_motor1.p1; xfit_motor1.p2; xfit_motor1.p3; xfit_motor1.p4; xfit_motor1.p5];
z_coeffs = [zfit_motor1.p1; zfit_motor1.p2; zfit_motor1.p3; zfit_motor1.p4; zfit_motor1.p5];

writematrix(fit_coeffs,'fitted_coeffs_v4.csv')

figure(3)
hold on
plot(xfit_motor1,output_fitv3_data(1,:),dl_fitv3,'.')

points = 1:res;
figure(5)
plot(points,output_fitv3_data(1,:),'.')
hold on; grid on;
plot(points,pred_fitv3_data(1,:),'.')

x_req = linspace(0,30,res);
z_pred = z_coeffs(1)*x_req.^4 + z_coeffs(2)*x_req.^3 + z_coeffs(3)*x_req.^2 + z_coeffs(4)*x_req.^1 + z_coeffs(5);

figure(6)
plot(zfit_motor1,output_fitv3_data(1,:),output_fitv3_data(3,:),'.')
grid on;
xlabel('x')
ylabel('z')


req = [x_req; zeros(1,res); z_pred];

error = norm(output_fitv3_data-req)

avg_error = error/res

writematrix(z_coeffs,'z_coeffs.csv')

figure(7)
plot3(output_fitv3_data(1,:),output_fitv3_data(2,:),output_fitv3_data(3,:),'.')
hold on; grid on; axis equal;
plot3(req(1,:),req(2,:),req(3,:),'.')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Data','Prediction')

