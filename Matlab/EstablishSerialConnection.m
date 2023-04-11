%Establish connection to Arduino Due and initiate 'serialdecObj'

clear all

a = arduino('/dev/tty.usbmodem11201','Due') %find arduino
serialdevObj = device(a,'SerialPort',1) %initiate serial port

%past possible lines for record keeping:
% arduinoObj = arduino("usbmodem11201","Due","Libraries",{'SPI','Serial','I2C'})
% serialdevObj = device(arduinoObj,'SerialPort',1)

% a = arduino('usbmodem1201','Due','Libraries','Serial','TraceOn',true)

% write(serialdevObj,[88 99 65]);
% numBytes = serialdevObj.NumBytesAvailable
% if numBytes > 0
%     response = read(serialdevObj,numBytes)
% end
