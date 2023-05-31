clear;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;

PROBE_ID_PENCOIL  = 0x3CC3F000;   % pen probe coil
PROBE_ID_SPINE    = 0x3D4C7400;   % Brook's spine coil

%read in coefficients and workspace
x_coeffs_tab = readtable("x_coeffs_full.csv");
z_coeffs_tab = readtable("zfits_poly3_full_circle.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs = table2array(z_coeffs_tab);
workspace_tab = readtable("known_surface.csv");
workspace = table2array(workspace_tab);

%read in aurora information
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);
spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);
% load tip file
penprobe_tip = load('tipcal.tip');

%open serial port
device = serialport("/dev/tty.usbmodem11301",115200);
flush(device);


    

%create point requests and conver to theta and r
% pen = pendata_in_modelspace(:,rep);
% dist = vecnorm(workspace-pen);
% [val,indx] = min(dist);
% req(:,rep) = workspace(:,indx)
% r_req = sqrt(req(1,rep)^2+req(2,rep)^2)
% 
% if r_req > 33
%     r_req = 33
% end
% 
% theta = atan2(req(2,rep),req(1,rep))
% if theta < 0
%     theta = 2*pi+theta
% end

theta = 3*pi/2;
r_req = 5;

%linearly interpolate between known theta curves
ax = interp1(x_coeffs(1,:),x_coeffs(2,:),theta,"linear");
bx = interp1(x_coeffs(1,:),x_coeffs(3,:),theta,"linear");
cx = interp1(x_coeffs(1,:),x_coeffs(4,:),theta,"linear");
dx = interp1(x_coeffs(1,:),x_coeffs(5,:),theta,"linear");
ex = interp1(x_coeffs(1,:),x_coeffs(6,:),theta,"linear");
az = interp1(z_coeffs(1,:),z_coeffs(2,:),theta,"linear");
bz = interp1(z_coeffs(1,:),z_coeffs(3,:),theta,"linear");
cz = interp1(z_coeffs(1,:),z_coeffs(4,:),theta,"linear");
dz = interp1(z_coeffs(1,:),z_coeffs(5,:),theta,"linear");

%create cable change array
dlplane = ax*r_req.^4 + bx*r_req.^3 + cx*r_req.^2 + dx*r_req.^1 + ex;
dl = [cos(theta)*dlplane; -cos(theta)*dlplane; sin(theta)*dlplane; -sin(theta)*dlplane];

%convert to PWM motor inputs
dtheta = dl./(disk_diameter/2);
dPWM = -dtheta*PWMrange/pi;
u(1:4,1) = setmid;
u(1:4) = u(1:4)+dPWM;

write(device,u(:),"uint16")
count = size(u(:),1);
response = read(device,count,"uint16")

