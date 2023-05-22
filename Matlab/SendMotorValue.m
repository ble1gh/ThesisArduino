%Brook Leigh
clear all

%Establish serial connection
device = serialport("/dev/tty.usbmodem1301",115200)
flush(device);

%Enter motor values
pos1 = 350;
% pos3 = 310;
% pos2 = 350;
% pos4 = 340;
motorvalue = [pos1 pos1 pos1 pos1]';
%motorvalue = [300 300 300 300]';
%motorvalue = [320 320 320 320]';


%Write to device and read response
write(device,motorvalue,"uint16")
count = size(motorvalue,1);

response = read(device,count,"uint16")'
%response = flip(response,2)'
assert(all((response == motorvalue)'),'Response does not match requested motor values. Motors are likely saturated')