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
dlplane = zeros(res_curve*iterations,res_theta);
zfits5 = zeros(6,res_theta);
xfits = zeros(5,res_theta);

%values to be used in calcs
mid = ones(res_curve*iterations,1)*setmid;
theta = [linspace(0,pi/2,10) linspace(pi/2,pi,10) linspace(pi,3*pi/2,10) linspace(3*pi/2,2*pi,10)];


figure(1)
hold on; grid on;
for i = 1:res_theta
    for k = 1:iterations
        points_record(:,:,k,i) = points_record_2D(:,(i-1)*res_curve+1:i*res_curve);
        x((k-1)*res_curve+1:k*res_curve,i) = points_record(1,:,k,i);
        y((k-1)*res_curve+1:k*res_curve,i) = points_record(2,:,k,i);
        z((k-1)*res_curve+1:k*res_curve,i) = points_record(3,:,k,i);
        r((k-1)*res_curve+1:k*res_curve,i) = sqrt(points_record(1,:,k,i).^2+points_record(2,:,k,i).^2);
        u1((k-1)*res_curve+1:k*res_curve,i) = points_record(4,:,k,i);
    end
    %find fits incorporating all 3 iterations
    % dl1(:,i) = -(u1(:,i)-mid)*pi/PWMrange*(disk_diameter/2);
    % dlplane(:,i) = dl1(:,i)./cos(theta(i));
    % [xfit_r, gofx] = fit(r(:,i),dlplane(:,i),'poly4');
    % xfits(:,i) = [xfit_r.p1; xfit_r.p2; xfit_r.p3; xfit_r.p4; xfit_r.p5];
    % plot(xfit_r,r(:,i),dlplane(:,i))
    % legend off
    [zfit_r, gofz] = fit(r(:,i),z(:,i),'poly5')
    zfits5(:,i) = [zfit_r.p1; zfit_r.p2; zfit_r.p3; zfit_r.p4; zfit_r.p5; zfit_r.p6];
    plot(zfit_r,r(:,i),z(:,i))
    legend off
end

zfits5 = [theta; zfits5];

zfits5_full_circle = [zfits5(:,1:9) zfits5(:,11:19) zfits5(:,21:29) zfits5(:,31:40)];

writematrix(zfits5_full_circle,'zfits_poly5_full_circle.csv')

%dlplane_onerep(1,:) = xfit_coeffs(1,i)*r_req.^4 + xfit_coeffs(2,i)*r_req.^3 + xfit_coeffs(3,i)*r_req.^2 + xfit_coeffs(4,i)*r_req.^1 + xfit_coeffs(5,i);

%writematrix(zfits,"zfits_3runs.csv")

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

%add the zero points to r for fitting. They already exist in the dl data
r0 = zeros(1,res_theta);
r = [r0; r(1:18,:); r0; r(19:18*2,:); r0; r(18*2+1:18*3,:)];

xfits5 = zeros(6,res_theta);

figure(2)
hold on; grid on;
for i = 1:res_theta
    [xfit_r, gofz] = fit(r(:,i),dl(:,i),'poly5');
    xfits5(:,i) = [xfit_r.p1; xfit_r.p2; xfit_r.p3; xfit_r.p4; xfit_r.p5; xfit_r.p6];
    plot(xfit_r,r(:,i),dl(:,i))
end
legend off

xfits5 = [theta; xfits5];

xfits5_full_circle = [xfits5(:,1:9) xfits5(:,11:19) xfits5(:,21:29) xfits5(:,31:40)];

writematrix(xfits5_full_circle,'xfits_poly5_full_circle.csv');

