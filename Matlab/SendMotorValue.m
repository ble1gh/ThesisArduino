%Brook Leigh
clear all

%Establish serial connection
device = serialport("/dev/tty.usbmodem11301",9600)
flush(device);

%Enter motor values
%motorvalue = [300 300 300 300];
motorvalue = 300;

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);
response = read(device,count(2),"uint16")