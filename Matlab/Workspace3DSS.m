%Brook Leigh
%Simulation of 3D reachable area for one segment and find motor inputs and
%send them to arduino
clear

L = 70e-3; % (m) from spine architecture 
d = 4e-3;
disk_diameter = 15e-3; % (m) diameter of the disk
PWMrange = 500-100;
res = 10;

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(3*L/(2*pi));
maxR = log10(2);
Rc = [logspace(minR,maxR,res) -logspace(maxR,minR,res)];

%2 alternative ways to get acceptable Rc values
% Rc = [linspace(3*L/(2*pi),2*L/pi,50) linspace(2*L/pi,.8,100) linspace(-.8,-2*L/pi,100) linspace(-2*L/pi,-3*L/(2*pi),50)]; 
% Rc = [linspace(3*L/(2*pi),2*L/pi,50) linspace(2*L/pi,.8,100)]; 

theta = linspace(0,pi,res);
phi = L./Rc;

%initiate state variable
x = zeros(3,length(theta)*length(Rc));
dl = zeros(4,length(theta)*length(Rc));

%for a given theta, an arc is defined along the range of possible bend
%radiuses which each give an (x,y,z) coordinate
l = length(Rc);
for i = 1:length(theta)
    x(1,((i-1)*l+1):i*l) = cos(theta(i))*Rc.*(1-cos(L./Rc));
    x(2,((i-1)*l+1):i*l) = sin(theta(i))*Rc.*(1-cos(L./Rc));
    x(3,((i-1)*l+1):i*l) = Rc.*sin(L./Rc);
    dl(1,((i-1)*l+1):i*l) = cos(theta(i)).*(phi.*(Rc-d)-L);
    dl(2,((i-1)*l+1):i*l) = cos(theta(i)).*(phi.*(Rc+d)-L);
    dl(3,((i-1)*l+1):i*l) = sin(theta(i)).*(phi.*(Rc-d)-L);
    dl(4,((i-1)*l+1):i*l) = sin(theta(i)).*(phi.*(Rc+d)-L);
end

%turn length change values into input values
u = ones(4,length(dl))*300;
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u = u+dPWM;

%plot result
figure(1)
plot3(x(1,:),x(2,:),x(3,:),'.')
grid on
title('Simulated Reachable Region Discretized, 3D')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-.07 .07 -.07 .07 0 .08])
%% 

%Establish serial connection
device = serialport("/dev/tty.usbmodem1301",115200)
flush(device);

%Enter motor values
motorvalue = [300 300 300 300];

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(2),"uint16")

%now do same for all u sets
for k = 1:length(u)
    motorvalue = u(:,k)
    k
    write(device,motorvalue,"uint16")
    count = size(motorvalue);
    response = read(device,count(1),"uint16")
    %w = waitforbuttonpress;
    pause(.5)
end

clear device