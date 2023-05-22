clear;

%physical system parameters
L = 55; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 340;

%resolution
res =18;

%Rc range
minR = log10(2*L/(pi));
maxR = log10(300);
Rc = [logspace(log10(3000),maxR,res/2) logspace(maxR,minR,res)];
phi = L./Rc;
theta = 0;

%initialize state and input vars
x = zeros(2,length(Rc));
dl = zeros(4,length(Rc)); 
u = ones(4,length(dl))*setmid;

%predicted xz
x(1,:) = -Rc.*(1-cos(L./Rc));
x(2,:) = Rc.*sin(L./Rc);

%cable length change
dl(1,:) = (phi.*(Rc-d*cos(theta))-L);
dl(2,:) = (phi.*(Rc+d*cos(theta))-L);
dl(3,:) = (phi.*(Rc-d*sin(theta))-L);
dl(4,:) = (phi.*(Rc+d*sin(theta))-L);

%convert to motor values
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u = u-dPWM;

figure(1)
plot(x(1,:),x(2,:),'.')
xlabel('x (mm)')
ylabel('z (mm)')
title('Expected Points')

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
writematrix(u,'modelb_input.csv')
writematrix(x,'modelb_predicted.csv')
writematrix(xstar,'modelb_output.csv')
writematrix(Rc,'Rc_data.csv')
writematrix(dl,'modelb_cable_dl.csv')

%disconnect serial ports
clear fser
clear device