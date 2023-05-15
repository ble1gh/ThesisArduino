%Brook Leigh
clear all

%Establish serial connection
device = serialport("/dev/tty.usbmodem11301",115200)
flush(device);

%Enter motor values
pos1 = 340;
% pos3 = 310;
% pos2 = 350;
% pos4 = 350;
motorvalue = [pos1 pos1 pos1 pos1]';
motorvalue = [300 300 300 300]';
%motorvalue = u(:,250)

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);

response = read(device,count(1),"uint16")