clear;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 340;

%resolution
res =18;

%Rc range
minR = log10(4*L/(pi));
maxR = log10(300);
Rc = [logspace(log10(3000),maxR,res/2) logspace(maxR,minR,res)];
phi = L./Rc;
theta = linspace(0,pi/2,10);

%initialize state and input vars
x = zeros(3,length(Rc)); 
dl = zeros(4,length(Rc)); 
u = ones(4,length(Rc))*setmid;

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

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
count = size(motorvalue);
response = read(device,count(1),"uint16")

%array to record actual positions
xstar = zeros(length(u),10,length(theta));

for i = 1:length(theta)
    %predicted xyz
    x(1,:) = cos(theta(i))*Rc.*(1-cos(L./Rc));
    x(2,:) = sin(theta(i))*Rc.*(1-cos(L./Rc));
    x(3,:) = Rc.*sin(L./Rc);
    
    %cable length change
    dl(1,:) = (phi.*(Rc-d*cos(theta(i)))-L);
    dl(2,:) = (phi.*(Rc+d*cos(theta(i)))-L);
    dl(3,:) = (phi.*(Rc-d*sin(theta(i)))-L);
    dl(4,:) = (phi.*(Rc+d*sin(theta(i)))-L);
    
    %convert to motor values
    dtheta = dl/(disk_diameter/2);
    dPWM = dtheta*PWMrange/pi;
    u = u-dPWM;
    
    figure(1)
    plot3(x(1,:),x(2,:),x(3,:),'.')
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('Expected Points')
    hold on; grid on; axis equal;
    patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);

    %now do same for all u sets
    for k = 1:length(u)
        motorvalue = round(u(:,k))
        k
        write(device,motorvalue,"uint16")
        count = size(motorvalue);
        response = read(device,count(1),"uint16")'
        %response = flip(response,2)'
        assert(all((response == motorvalue)'),'Response does not match requested motor values. Motors are likely saturated')
        pause(2)
    
        pkt = [];
        while(isempty(pkt))
            pkt = getAuroraPacket(fser,0.1);
        end
        xstar(k,:,i) = pkt;
        i
    end
    
    %return to neutral
    motorvalue = [setmid setmid setmid setmid]';
    write(device,motorvalue,"uint16")
    count = size(motorvalue);
    response = read(device,count(1),"uint16")

    %write data to csv files
    input_file = ['input_theta_' num2str(theta(i),'%.2f') '.csv'];
    writematrix(u,input_file);
    predicted_file = ['predicted_theta_' num2str(theta(i),'%.2f') '.csv'];
    writematrix(x,predicted_file);
    output_file = ['output_theta_' num2str(theta(i),'%.2f') '.csv'];
    writematrix(xstar(:,:,i),output_file);
end


%disconnect serial ports
clear fser
clear device