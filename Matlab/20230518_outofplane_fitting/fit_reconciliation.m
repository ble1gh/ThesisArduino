clear; close all;

%load in fit coefficients
xfit_tab = readtable("xfit_coeffs.csv");
xfit_coeffs = table2array(xfit_tab);
zfit_tab = readtable("zfit_coeffs.csv");
zfit_coeffs = table2array(zfit_tab);


%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;

%set frequency of sampling
res_curve = 18;
res_theta = 10;
iterations = 5;

r_req = linspace(0,30,res_curve);
dlplane = zeros(iterations,res_curve,res_theta);
theta = linspace(0,pi/2,res_theta);
final_coeffs = ones(5,res_theta);

for i = 1:res_theta
    for j = 2:iterations+1
        dlplane(j-1,:,i) = xfit_coeffs(1,(i-1)*(iterations+1)+j)*r_req.^4 + xfit_coeffs(2,(i-1)*(iterations+1)+j)*r_req.^3 + xfit_coeffs(3,(i-1)*(iterations+1)+j)*r_req.^2 + xfit_coeffs(4,(i-1)*(iterations+1)+j)*r_req.^1 + xfit_coeffs(5,(i-1)*(iterations+1)+j);
    end
    final_coeffs(:,i) = [xfit_coeffs(1,i*6); xfit_coeffs(2,i*6); xfit_coeffs(3,i*6); xfit_coeffs(4,i*6); xfit_coeffs(5,i*6)];
end

for i = 1:res_theta-3
    for j = 1:iterations
        figure(i)
        hold on; grid on;
        plot(r_req,dlplane(j,:,i))
        xlabel('r (mm)')
        ylabel('dl In Plane (mm)')
        legend('1','2','3','4','5')
    end
end

% for i = 1:res_theta-1
% 
% end