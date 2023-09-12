%Brook Leigh
%Mapping of 3D reachable area for one segment. Constant Curvature method
%find motor inputs and send them to arduino
clear; close all;

L = 56; % (mm) from spine architecture 
d = 4; % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
res = 25;
setmid = 350;

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(2.5*L/(pi));
maxR = log10(200);
Rc = [logspace(minR,maxR,res/2)];

theta = linspace(0,2*pi,res);
phi = L./Rc;

%initiate state variable
x = zeros(3,length(theta)*length(Rc));
dl = zeros(4,length(theta)*length(Rc));

% for a given Rc, an arc is defined along the range of possible
% thetas which each give an (x,y,z) coordinate, and thus a ring
l = length(Rc);
for i = 1:length(theta)
    x(1,((i-1)*l+1):i*l) = cos(theta(i))*Rc.*(1-cos(L./Rc));
    x(2,((i-1)*l+1):i*l) = sin(theta(i))*Rc.*(1-cos(L./Rc));
    x(3,((i-1)*l+1):i*l) = Rc.*sin(L./Rc);

    dl(1,((i-1)*l+1):i*l) = (phi.*(Rc-d*cos(theta(i)))-L);
    dl(2,((i-1)*l+1):i*l) = (phi.*(Rc+d*cos(theta(i)))-L);
    dl(3,((i-1)*l+1):i*l) = (phi.*(Rc-d*sin(theta(i)))-L);
    dl(4,((i-1)*l+1):i*l) = (phi.*(Rc+d*sin(theta(i)))-L);
end

%turn length change values into input values
u = ones(4,length(dl))*setmid;
dtheta = dl/(disk_diameter/2);
dPWM = -dtheta*PWMrange/pi;
u = u+dPWM;

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

%plot result
figure(1)
plot3(x(1,:),x(2,:),x(3,:),'.')
grid on; hold on;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
title('Predicted Path')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

%visualization of motor values
% figure(2)
% t = 1:length(u);
% plot(theta,u(1,:),theta,u(2,:),theta,u(3,:),theta,u(4,:))
% title('Motor Values')
% legend('1','2','3','4')
% xlabel('theta')
% ylabel('PWM Value')
% grid on

% figure(3)
% t = 1:length(u);
% plot(theta,dl(1,:),theta,dl(2,:),theta,dl(3,:),theta,dl(4,:),theta,(phi*(Rc-d*cos(theta))-L))
% legend('1','2','3','4','(phi*(Rc-d*cos(theta))-L)')
% xlabel('theta')
% ylabel('dl')
% grid on

runtime_minutes = length(u)*2.5/60
%%

%read in aurora information
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);
spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);

PROBE_ID_PENCOIL  = 0x3CC3F000;   % pen probe coil
PROBE_ID_SPINE    = 0x3D4C7400;   % Brook's spine coil

%Establish serial connection to Arduino
device = serialport("/dev/tty.usbmodem11301",115200)
flush(device);

%Establish serial connection to Aurora
fser = serialport('/dev/tty.usbserial-1120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);

%Enter motor values
motorvalue = [setmid setmid setmid setmid]';

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue,1);
response = read(device,count,"uint16")

%array to record actual positions
spinedata_in_modelspace = zeros(3,length(u));

%now do same for all u sets
for j = 1:length(u)
    motorvalue = u(:,j)

    write(device,motorvalue,"uint16")
    count = size(motorvalue,1);
    response = read(device,count,"uint16")
    pause(1)

    pkt = [];
    while(isempty(pkt))
        requestAuroraPacket(fser,[PROBE_ID_SPINE]);
        pkt = getAuroraPacket(fser,0.2);
    end

    assert(mod(length(pkt)-1,9) == 0,'Incorrect packet size!');
    num_tforms = (length(pkt)-1)/9;
    all_tfs = [];
    for tf_idx = 1:num_tforms
        all_tfs(tf_idx).sn = pkt(2+9*(tf_idx-1));
        all_tfs(tf_idx).T_coil_to_aurora = eye(4);
        all_tfs(tf_idx).T_coil_to_aurora(1:3,1:3) = quat2matrix(pkt((3:6)+9*(tf_idx-1)));
        all_tfs(tf_idx).T_coil_to_aurora(1:3,4) = pkt((7:9)+9*(tf_idx-1));
        all_tfs(tf_idx).error = pkt(10+9*(tf_idx-1));
    end
    if(num_tforms)
        all_tfs.T_coil_to_aurora
    else
        warning('No transforms returned!');
    end

    spinedata_in_auroraspace = zeros(3,1);

    for i = 1:length(all_tfs)
        if all_tfs(i).sn == PROBE_ID_SPINE
            spinedata_in_auroraspace = (all_tfs(i).T_coil_to_aurora(1:3,4) + all_tfs(i).T_coil_to_aurora(1:3,1:3)*spinetip_in_coilspace);
        end
    end
    
    spinedata_in_modelspace(:,j) = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);
end

%return to neutral
motorvalue = [setmid setmid setmid setmid]';
write(device,motorvalue,"uint16")
count = size(motorvalue,1);
response = read(device,count,"uint16")

points_record_test = zeros(7,length(u));
points_record_test(1:3,:) = spinedata_in_modelspace;
points_record_test(4:7,:) = u;
writematrix(points_record_test,'points_record_model_test.csv');
writematrix(x,'points_request_model_test.csv');

figure(1)
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','Color',"#77AC30");
legend('Requests','Position Data')
hold on; grid on; axis equal;
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

clear device