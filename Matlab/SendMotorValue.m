%Brook Leigh
clear all

%Establish serial connection
device = serialport("/dev/tty.usbmodem11301",115200)
flush(device);

%Enter motor values
pos1 = 300;
pos3 = 300;
pos2 = 300 + (300-pos1);
pos4 = 300 + (300-pos3);
motorvalue = [pos1 pos2 pos3 pos4]';
% motorvalue = [300 300 300 300]';
%motorvalue = u(:,250)

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);

response = read(device,count(1),"uint16")