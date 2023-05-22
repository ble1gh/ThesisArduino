clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;
theta = 0;
res = 18;

fit_coeffs_tab = readtable("fitted_coeffs_v3.csv");
fit_coeffs = table2array(fit_coeffs_tab);

k =2;

x_req = linspace(0,30,res);
dl1 = fit_coeffs(1,k)*x_req.^4 + fit_coeffs(2,k)*x_req.^3 + fit_coeffs(3,k)*x_req.^2 + fit_coeffs(4,k)*x_req.^1 + fit_coeffs(5,k);
dtheta1 = dl1/(disk_diameter/2);
dPWM1 = -dtheta1*PWMrange/pi;
dPWM2 = -dPWM1;

u = ones(4,length(x_req))*setmid;
u(1,:) = u(1,:)+dPWM1;
u(2,:) = u(2,:)+dPWM2;

figure(1)
plot(x_req,u(1,:),'.')
hold on; grid on;
plot(x_req,u(2,:),'.')
plot(x_req,u(3,:),'.')
plot(x_req,u(4,:),'.')
%%

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
xstar = zeros(length(u),10);

%now do same for all u sets
for k = 1:length(u)
    motorvalue = u(:,k)
    k
    write(device,motorvalue,"uint16")
    count = size(motorvalue,1);
    response = read(device,count,"uint16")
    pause(2.5)

    pkt = [];
    while(isempty(pkt))
        pkt = getAuroraPacket(fser,0.1);
    end
    xstar(k,:) = pkt;
end

%return to neutral
motorvalue = [setmid setmid setmid setmid]';
write(device,motorvalue,"uint16")
count = size(motorvalue,1);
response = read(device,count,"uint16")

%write data to csv files
writematrix(u,'fitted_v3_input.csv')
writematrix(x_req,'fitted_v3_predicted.csv')
writematrix(xstar,'fitted_v3_output.csv')
writematrix(dl1,'fitted_v3_cable_dl.csv')

%disconnect serial ports
clear fser
clear device