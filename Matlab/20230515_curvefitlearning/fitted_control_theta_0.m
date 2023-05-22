clear; close all;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 340;
theta = 0;

%read in data
Rc_tab = readtable('Rc_range.csv');
Rc_data = -table2array(Rc_tab);
predicted_tab = readtable('predicted_theta_0.00.csv');
predicted_data = table2array(predicted_tab);
output_tab = readtable('output_theta_0.00.csv');
output_data = table2array(output_tab);
pos_data = output_data(:,7:9)';

%calculate error
model_error = norm(pos_data-predicted_data)

%curve fit
[xfit_Rc, gofx] = fit(pos_data(1,:)',Rc_data','exp2');

figure(1)
plot(xfit_Rc,pos_data(1,:),Rc_data,'.')
hold on; grid on;
plot(predicted_data(1,:),Rc_data)
legend('Data','Curve Fit','Prediction','Location','east')
xlabel('x distance from center (mm)')
ylabel('Rc')
title('Rc vs x, Prediction vs Reality')

%generate new points to test
R_fit = [xfit_Rc.a; xfit_Rc.b; xfit_Rc.c; xfit_Rc.d];
x = linspace(3,20,20);
R = R_fit(1)*exp(R_fit(2)*x) + R_fit(3)*exp(R_fit(4)*x);

%initialize prediction and input vars
pred = zeros(2,length(R));
dl = zeros(4,length(R)); 
u = ones(4,length(dl))*setmid;

%predicted xz
pred(1,:) = R.*(1-cos(L./R));
pred(2,:) = R.*sin(L./R);

%cable length change
dl(1,:) = (L./R.*(R-d*cos(theta))-L);
dl(2,:) = (L./R.*(R+d*cos(theta))-L);
dl(3,:) = (L./R.*(R-d*sin(theta))-L);
dl(4,:) = (L./R.*(R+d*sin(theta))-L);

%convert to motor values
dtheta = dl/(disk_diameter/2);
dPWM = dtheta*PWMrange/pi;
u = u-dPWM;

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
count = size(motorvalue);
response = read(device,count(1),"uint16")

%write data to csv files
writematrix(u,'fitted_input.csv')
writematrix(x,'fitted_predicted.csv')
writematrix(xstar,'fitted_output.csv')

%disconnect serial ports
clear fser
clear device
