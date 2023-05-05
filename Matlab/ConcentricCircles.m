%Brook Leigh
%Simulation of 3D reachable area for one segment and find motor inputs and
%send them to arduino
clear

L = 70; % (mm) from spine architecture 
d = 4; % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
res = 25;

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(2*L/(pi));
maxR = log10(150);
Rc = [logspace(minR,maxR,res/2)];

theta = linspace(0,2*pi,res);
phi = L./Rc;

% %One Circle
% Rc = Rc(1);
% phi = phi(1);

% %initiate state and output vars
% x = zeros(3,length(theta));
% dl = zeros(4,length(theta));
% 
% x(1,:) = cos(theta)*Rc*(1-cos(L/Rc));
% x(2,:) = sin(theta)*Rc*(1-cos(L/Rc));
% x(3,:) = Rc*sin(L/Rc);
% 
% dl(1,:) = (phi*(Rc-d*cos(theta))-L);
% dl(2,:) = (phi*(Rc+d*cos(theta))-L);
% dl(3,:) = (phi*(Rc-d*sin(theta))-L);
% dl(4,:) = (phi*(Rc+d*sin(theta))-L);

%initiate state variable
x = zeros(3,length(theta)*length(Rc));
dl = zeros(4,length(theta)*length(Rc));

% for a given Rc, an arc is defined along the range of possible
% thetas which each give an (x,y,z) coordinate, and thus a ring
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


%plot result
figure(1)
plot3(x(1,:),x(2,:),x(3,:))
grid on
title('Predicted Path')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
% axis([-.07 .07 -.07 .07 0 .08])

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


%% 

%Establish serial connection
device = serialport("/dev/tty.usbmodem11301",115200)
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
    pause(2.5)
end

motorvalue = [300 300 300 300]';

write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(1),"uint16")

clear device