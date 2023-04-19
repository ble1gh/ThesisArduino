%Brook Leigh
clear all

%Establish serial connection
device = serialport("/dev/tty.usbmodem11301",115200)
flush(device);

%Enter motor values
motorvalue = [300 299 298 297];
%motorvalue = 300;

%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue);

response = read(device,count(2),"uint16")