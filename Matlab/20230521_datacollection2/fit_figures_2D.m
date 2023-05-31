clear; close all;

%load full circle data history
points_record_tab = readtable('points_record_curvefit_test.csv');
points_record_2D = table2array(points_record_tab);
x_coeffs_tab = readtable('x_coeffs_full.csv');
x_coeffs = table2array(x_coeffs_tab);

dl_all_iterations_Q1_tab = readtable("dl_all_iterations_Q1.csv");
dl_all_iterations_Q1 = table2array(dl_all_iterations_Q1_tab);
dl_all_iterations_Q2_tab = readtable("dl_all_iterations_Q2.csv");
dl_all_iterations_Q2 = table2array(dl_all_iterations_Q2_tab);
dl_all_iterations_Q3_tab = readtable("dl_all_iterations_Q3.csv");
dl_all_iterations_Q3 = table2array(dl_all_iterations_Q3_tab);
dl_all_iterations_Q4_tab = readtable("dl_all_iterations_Q4.csv");
dl_all_iterations_Q4 = table2array(dl_all_iterations_Q4_tab);

dl2D = [dl_all_iterations_Q1 dl_all_iterations_Q2 dl_all_iterations_Q3 dl_all_iterations_Q4];

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

%initiate vars
points_record = zeros(7,res_curve,iterations,res_theta);
r = zeros(res_curve*iterations,res_theta);
x = zeros(res_curve*iterations,res_theta);
y = zeros(res_curve*iterations,res_theta);
z = zeros(res_curve*iterations,res_theta);
u1 = zeros(res_curve*iterations,res_theta);
dl1 = zeros(res_curve*iterations,res_theta);
dlplane = zeros(res_curve,iterations,res_theta);
zfits5 = zeros(6,res_theta);
xfits = zeros(5,res_theta);

%values to be used in calcs
mid = ones(res_curve*iterations,1)*setmid;
theta = [linspace(0,pi/2,10) linspace(pi/2,pi,10) linspace(pi,3*pi/2,10) linspace(3*pi/2,2*pi,10)];

dl = zeros(57,res_theta);

for i = 1:10
    dl(1:57,i) = dl_all_iterations_Q1(57*(i-1)+1:i*57)';
end

for i = 1:10
    dl(1:57,10+i) = dl_all_iterations_Q2(57*(i-1)+1:i*57)';
end

for i = 1:10
    dl(1:57,20+i) = dl_all_iterations_Q3(57*(i-1)+1:i*57)';
end

for i = 1:10
    dl(1:57,30+i) = dl_all_iterations_Q4(57*(i-1)+1:i*57)';
end

for i = 1:res_theta
    for k = 1:iterations
        points_record(:,:,k,i) = points_record_2D(:,(i-1)*res_curve+1:i*res_curve);
        dlplane(:,k,i) = dl2D(:,(i-1)*(res_curve+1)+2:i*(res_curve+1))';
        x((k-1)*res_curve+1:k*res_curve,i) = points_record(1,:,k,i);
        y((k-1)*res_curve+1:k*res_curve,i) = points_record(2,:,k,i);
        z((k-1)*res_curve+1:k*res_curve,i) = points_record(3,:,k,i);
        r((k-1)*res_curve+1:k*res_curve,i) = sqrt(points_record(1,:,k,i).^2+points_record(2,:,k,i).^2);
        u1((k-1)*res_curve+1:k*res_curve,i) = points_record(4,:,k,i);
    end
end

[xfit, gofx] = fit(points_record(1,:,1,1)',dlplane(:,1,1),'poly4')
figure(2)
plot(xfit,points_record(1,:,1,1),dlplane(:,1,1),'.')
grid on;
xlabel('r (mm)')
ylabel('Cable Length Change (mm)')


[xfit, gofx] = fit(r(18*2+1:18*3,33),dlplane(:,3,33),'poly4')

figure(1)
plot(xfit,r(18*2+1:18*3,33),dlplane(:,3,33),'.')
grid on;
xlabel('r (mm)')
ylabel('Cable Length Change (mm)')


%add the zero points to r for fitting. They already exist in the dl data
r0 = zeros(1,res_theta);
%r = [r0; r(1:18,:); r0; r(19:18*2,:); r0; r(18*2+1:18*3,:)];

%These values give a turn angle (theta) and a bend amount (Rc)
res = 500;
minR = log10(2.5*L/(pi));
maxR = log10(10000);
Rc = [logspace(minR,maxR,res)];
phi = L./Rc;

rmodel = Rc.*(1-cos(L./Rc));
zmodel = Rc.*sin(L./Rc);
dlmodel = (phi.*(Rc-d)-L);

%other interesting spot is theta(11) = pi/2

figure(1)
hold on;
plot(rmodel,dlmodel,'Color','k','LineWidth',2)
title('Cable Length Change vs Radial Displacement, theta = 5.3','FontSize',18)
xlabel('r (mm)')
ylabel('dl (mm)')
legend('data','4th order polynomial fit','model')

figure(2)
hold on;
plot(rmodel,dlmodel,'k','LineWidth',2)
title('Cable Length Change vs Radial Displacement, xz Plane','FontSize',18)
xlabel('r (mm)')
ylabel('dl (mm)')
legend('data','4th order polynomial fit','model')


[zfit, gofz] = fit(points_record(1,:,1,1)',points_record(3,:,1,1)','poly4')

figure(3)
plot(zfit,points_record(1,:,1,1)',points_record(3,:,1,1)')
hold on; grid on;
plot(rmodel,zmodel,'k','LineWidth',2)
title('Height vs Radial Displacement, xz Plane','FontSize',18)
xlabel('r (mm)')
ylabel('z (mm)')
legend('data','3rd order polynomial fit','model')
axis([0 35 0 60])

[xfit, gofx] = fit(r(18*2+1:18*3,33),z(18*2+1:18*3,33),'poly4')

figure(4)
plot(zfit,r(18*2+1:18*3,33),z(18*2+1:18*3,33),'.')
title('Cable Length Change vs Radial Displacement, theta = 5.3','FontSize',18)
hold on; grid on;
plot(rmodel,zmodel,'k','LineWidth',2)
xlabel('r (mm)')
ylabel('z (mm)')
legend('data','3rd order polynomial fit','model')

figure(5)
plot(rmodel,zmodel,'LineWidth',2)
grid on;
xlabel('r (mm)')
ylabel('z (mm)')
title('Model Prediction, r vs z In Plane','FontSize',18)
axis([0 35 0 60])

figure(6)
plot(rmodel,dlmodel,'LineWidth',2)
grid on;
xlabel('r (mm)')
ylabel('dl (mm)')
title('Model Prediction, r vs Cable Length Change','FontSize',18)
