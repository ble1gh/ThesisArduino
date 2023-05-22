clear;

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

%read in coefficients
x_coeffs_tab = readtable("x_coeffs_full.csv");
z_coeffs_tab = readtable("z_coeffs_full.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs = table2array(z_coeffs_tab);

%read in points array for position requesting
points_record_tab = readtable("points_record_circle.csv");
points_record = table2array(points_record_tab);

%create point requests
%establish requests, and resulting cable change
theta = linspace(0,2*pi,res_theta);
r_req = linspace(0,28,res_curve);

%linearly interpolate between known theta curves
ax = interp1(x_coeffs(1,:),x_coeffs(2,:),theta,"linear");
bx = interp1(x_coeffs(1,:),x_coeffs(3,:),theta,"linear");
cx = interp1(x_coeffs(1,:),x_coeffs(4,:),theta,"linear");
dx = interp1(x_coeffs(1,:),x_coeffs(5,:),theta,"linear");
ex = interp1(x_coeffs(1,:),x_coeffs(6,:),theta,"linear");
% az = interp1(z_coeffs(1,:),z_coeffs(2,:),theta,"linear");
% bz = interp1(z_coeffs(1,:),z_coeffs(3,:),theta,"linear");
% cz = interp1(z_coeffs(1,:),z_coeffs(4,:),theta,"linear");
% dz = interp1(z_coeffs(1,:),z_coeffs(5,:),theta,"linear");
% ez = interp1(z_coeffs(1,:),z_coeffs(6,:),theta,"linear");


% pos_req = zeros(3,res_curve*res_theta);
% for i = 1:res_theta
%     pos_req(:,(i-1)*res_curve+1:i*res_curve) = [r_req.*cos(theta(i)); r_req.*sin(theta(i)); 
%     az(i)*r_req.^4 + bz(i)*r_req.^3 + cz(i)*r_req.^2 + dz(i)*r_req.^1 + ez(i)];
% end


%create 2D array of pos_req for ploting
query = zeros(2,res_curve*res_theta);
for i = 1:res_theta
    query(:,(i-1)*res_curve+1:i*res_curve) = [r_req.*cos(theta(i)); r_req.*sin(theta(i))];
end

%use past data to interpolate z value of x,y pairs
x = sortrows(points_record(1,:)')';
y = sortrows(points_record(2,:)')';
[X,Y] = meshgrid(x,y);
Z = griddata(points_record(1,:),points_record(2,:),points_record(3,:),X,Y);
zq = interp2(X,Y,Z,query(1,:),query(2,:),"linear");
query = [query(1,:); query(2,:); zq];

figure(1)
hold on; grid on; axis equal;
plot3(query(1,:),query(2,:),query(3,:),'.')

%allocate variables for length change
dlplane = zeros(1,res_curve*res_theta);
dl = zeros(4,res_curve*res_theta);

%create cable change array
for i = 1:res_theta
    dlplane((i-1)*res_curve+1:i*res_curve) = ax(i)*r_req.^4 + bx(i)*r_req.^3 + cx(i)*r_req.^2 + dx(i)*r_req.^1 + ex(i);
    dl(1,(i-1)*res_curve+1:i*res_curve) = cos(theta(i))*dlplane((i-1)*res_curve+1:i*res_curve);
    dl(2,(i-1)*res_curve+1:i*res_curve) = -cos(theta(i))*dlplane((i-1)*res_curve+1:i*res_curve);
    dl(3,(i-1)*res_curve+1:i*res_curve) = sin(theta(i))*dlplane((i-1)*res_curve+1:i*res_curve);
    dl(4,(i-1)*res_curve+1:i*res_curve) = -sin(theta(i))*dlplane((i-1)*res_curve+1:i*res_curve);
end

%convert to PWM motor inputs
dtheta = dl./(disk_diameter/2);
dPWM = -dtheta*PWMrange/pi;
u = ones(4,res_curve*res_theta)*setmid;
u = u+dPWM;

%read in aurora information
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);
spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);
%spine tip in model space
spine_tip_modelspace = [0,0,L]';

%%
%Establish serial connection to Arduino
device = serialport("/dev/tty.usbmodem1301",115200)
flush(device);

%Establish serial connection to Aurora
fser = serialport('/dev/tty.usbserial-120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);

%Enter motor values
motorvalue = [setmid setmid setmid setmid]';

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue,1);
response = read(device,count,"uint16")

%array to record actual positions
xstar = zeros(length(u),10);

%now do same for all u sets
for j = 1:length(u)
    motorvalue = u(:,j)

    write(device,motorvalue,"uint16")
    count = size(motorvalue,1);
    response = read(device,count,"uint16")
    pause(1)

    pkt = [];
    while(isempty(pkt))
        pkt = getAuroraPacket(fser,0.1);
    end
    xstar(j,:) = pkt;
end

%return to neutral
motorvalue = [setmid setmid setmid setmid]';
write(device,motorvalue,"uint16")
count = size(motorvalue,1);
response = read(device,count,"uint16")

for coil_idx = 1:size(xstar,1)
    q = xstar(coil_idx,3:6);
    t = xstar(coil_idx,7:9);
    spinedata_in_auroraspace(:,coil_idx) = (t + quatrotate(q,spinetip_in_coilspace'))';
end

spinedata_in_modelspace = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);
points_record_test = zeros(7,res_theta*res_curve);
points_record_test(1:3,:) = spinedata_in_modelspace;
points_record_test(4:7,:) = u;
writematrix(points_record_test,'points_record_curvefit_test.csv');
writematrix(query,'points_request_curvefit_test.csv');

figure(1)
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','Color',"#77AC30");
legend('Requests','Position Data')
hold on; grid on; axis equal;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

