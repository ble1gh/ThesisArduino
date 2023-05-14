%Brook Leigh
%Comes up with a list of (x,y,z) points and corresponding motor values,
%then sends motor values to arduino, records actual position, and compares
%to predicted/requested
clear

L = 56; % (mm) from spine architecture 
d = 4; % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
res = 10;

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(2*L/(pi));
maxR = log10(150);
Rc = [logspace(minR,maxR,res/2)];

theta = linspace(0,2*pi,res);
phi = L./Rc;

%All the circles

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
u_model = ones(4,length(dl))*300;
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u_model = u_model+dPWM;


%plot result
figure(2)
plot3(x(1,:),x(2,:),x(3,:),'.','color','b')
grid on
title('Predicted Path')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-70 70 -70 70 0 80])

%load data into var
history_tab = readtable('output_input_history_0427.csv','Delimiter',',');
ustar = table2array(history_tab);

xstar = sortrows(ustar(1,:)')';
ystar = sortrows(ustar(2,:)')';
[X,Y] = meshgrid(xstar,ystar);
Z = griddata(ustar(1,:),ustar(2,:),ustar(3,:),X,Y);

xq = zeros(1,length(x));
yq = zeros(1,length(x));
zq = zeros(1,length(x));

u_interp = zeros(4,length(x));

for i = 1:length(x)
    % here's a query point (used an existing point to test)
    xq(i) = x(1,i);
    yq(i) = x(2,i);
    
    % find associated z value for (x,y) pair
    
    zq(i) = interp2(X,Y,Z,xq(i),yq(i),"linear");
    
    %create interpolant objects for each motor given (x,y,z) position
    F1 = scatteredInterpolant(ustar(1:3,:)',ustar(4,:)','linear','linear');
    F2 = scatteredInterpolant(ustar(1:3,:)',ustar(5,:)','linear','linear');
    F3 = scatteredInterpolant(ustar(1:3,:)',ustar(6,:)','linear','linear');
    F4 = scatteredInterpolant(ustar(1:3,:)',ustar(7,:)','linear','linear');
    
    %evaluate the interpolant objects to find motor values given a (x,y,z) pos
    % m1 = round(F1(xq,yq,zq))
    % m2 = round(F2(xq,yq,zq))
    % m3 = round(F3(xq,yq,zq))
    % m4 = round(F4(xq,yq,zq))
    
    m1 = F1(xq(i),yq(i),zq(i));
    m2 = F2(xq(i),yq(i),zq(i));
    m3 = F3(xq(i),yq(i),zq(i));
    m4 = F4(xq(i),yq(i),zq(i));

    u_interp(:,i) = [m1 m2 m3 m4]';
end

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
xstar = zeros(length(u_model),10);

%now do same for all u sets
for k = 1:length(u_interp)
    motorvalue = u_interp(:,k)
    k
    write(device,motorvalue,"uint16")
    count = size(motorvalue);
    response = read(device,count(1),"uint16")
    pause(2.5)

    pkt = [];
    while(isempty(pkt))
        pkt = getAuroraPacket(fser,0.1);
    end
    xstar(k,:) = pkt;
end

%return to neutral
motorvalue = [300 300 300 300]';
write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(1),"uint16")

% 
%write data to csv files
writematrix(u_model,'input_interp.csv')
writematrix(x,'predicted_interp.csv')
writematrix(xstar,'output_interp.csv')

%disconnect serial ports
clear fser
clear device