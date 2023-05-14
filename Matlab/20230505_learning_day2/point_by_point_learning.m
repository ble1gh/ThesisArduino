% Brook Leigh
% going each time from a centered position, it goes to a point, then uses 
% real-time position data to narrow in on a point until it hits it

clear

L = 56; % (mm) from spine architecture 
d = 4; % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
res = 25;

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(2*L/(pi));
maxR = log10(150);
Rc = [logspace(minR,maxR,res/2)];

theta = linspace(2*pi,0,res);
phi = L./Rc;

%initiate state variable
x = zeros(3,length(theta)*length(Rc));
dl = zeros(4,length(theta)*length(Rc));

%for a given theta, an arc is defined along the range of possible bend
%radiuses which each give an (x,y,z) coordinate
l = length(theta);
for i = 1:length(Rc)
    x(1,((i-1)*l+1):i*l) = cos(theta)*Rc(i)*(1-cos(L/Rc(i)));
    x(2,((i-1)*l+1):i*l) = sin(theta)*Rc(i)*(1-cos(L/Rc(i)));
    x(3,((i-1)*l+1):i*l) = Rc(i)*sin(L/Rc(i));

    dl(1,((i-1)*l+1):i*l) = (phi(i)*(Rc(i)-d*cos(theta))-L);
    dl(2,((i-1)*l+1):i*l) = (phi(i)*(Rc(i)+d*cos(theta))-L);
    dl(3,((i-1)*l+1):i*l) = (phi(i)*(Rc(i)-d*sin(theta))-L);
    dl(4,((i-1)*l+1):i*l) = (phi(i)*(Rc(i)+d*sin(theta))-L);
end

%turn length change values into input values
u = ones(4,length(dl))*300;
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u = u+dPWM;

bent_motorvals = [x; u];
eq_pos = ones(size(bent_motorvals));
eq_pos(1:3,:) = eq_pos(1:3,:).*[0 0 L]';
eq_pos(4:7,:) = eq_pos(4:7,:).*[300 300 300 300]';
request_motorvals = [bent_motorvals zeros(size(eq_pos))];

for k = 1:length(bent_motorvals)
    request_motorvals(:,k*2-1) = eq_pos(:,k);
    request_motorvals(:,k*2) = bent_motorvals(:,k);
end

%plot result
figure(1)
plot3(x(1,:),x(2,:),x(3,:),'.','color','b')
grid on
title('Predicted Path')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-70 70 -70 70 0 80])

runtime_minutes = length(request_motorvals)*2.5/60

%%
%get current registration transformation matrices
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv','Delimiter',',');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);

spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv','Delimiter',',');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);

%% 

%Establish serial connection to Arduino
device = serialport("/dev/tty.usbmodem11301",115200)
flush(device);

%Establish serial connection to Aurora
fser = serialport('/dev/tty.usbserial-1120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);

%Enter motor values
motorvalue = [300 300 300 300];

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(2),"uint16")

%array to record actual positions
xstar_raw = zeros(length(request_motorvals),10);
xstar_auroraspace = zeros(3,length(request_motorvals));
xstar = zeros(3,length(request_motorvals));

%now do same for all u sets
for k = 1:length(request_motorvals)
    %send requested position
    motorvalue = request_motorvals(4:7,k)
    k
    write(device,motorvalue,"uint16")
    count = size(motorvalue);
    response = read(device,count(1),"uint16")
    pause(2.5)
    
    %request packet from aurora
    pkt = [];
    while(isempty(pkt))
        pkt = getAuroraPacket(fser,0.1);
    end

    %record raw data and transform to modelspace
    xstar_raw(k,:) = pkt;
    q = pkt(3:6);
    t = pkt(7:9);
    xstar_auroraspace(:,k) = (t + quatrotate(q,spinetip_in_coilspace'))';
    xstar(:,k) = hTF(xstar_auroraspace,TF_aurora_to_model,0);
    
    

end

%return to neutral
motorvalue = [300 300 300 300]';
write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(1),"uint16")


%disconnect serial ports
clear fser
clear device