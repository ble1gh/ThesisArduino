clear; close all;

%load full circle data history
points_record_tab = readtable('points_record_circle.csv');
points_record_2D = table2array(points_record_tab);
x_coeffs_tab = readtable('x_coeffs_full.csv');
x_coeffs = table2array(x_coeffs_tab);

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
zfits = zeros(5,res_theta);
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
    [zfit_r, gofz] = fit(r(:,i),z(:,i),'poly4');
    zfits(:,i) = [zfit_r.p1; zfit_r.p2; zfit_r.p3; zfit_r.p4; zfit_r.p5];
end

%dlplane_onerep(1,:) = xfit_coeffs(1,i)*r_req.^4 + xfit_coeffs(2,i)*r_req.^3 + xfit_coeffs(3,i)*r_req.^2 + xfit_coeffs(4,i)*r_req.^1 + xfit_coeffs(5,i);

zfits = [theta; zfits];

writematrix(zfits,"zfits_3runs.csv")