%recieve a packet from a constant stream of position data from Aurora
%through pc linux
fser = serialport('/dev/tty.usbserial-1120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);

pkt = [];
while(isempty(pkt))
    pkt = getAuroraPacket(fser,0.1);
end
%xstar(k,:) = pkt;

pkt

clear fser