%Brook Leigh
%Simulation of 2D reachable area for one segment
clear

L = 70e-3; % (m) from spine architecture 

Rc = [linspace(3*L/(2*pi),2*L/pi,100) linspace(2*L/pi,1,100) linspace(-1,-2*L/pi,100) linspace(-2*L/pi,-3*L/(2*pi),100)]; %just finding a reasonable range of values based on how far I want the spine to bend
%Rc = [linspace(3*L/(2*pi),2*L/pi,100) linspace(2*L/pi,1,100)];

x = -Rc.*(1-cos(L./Rc));
z = Rc.*sin(L./Rc);

% xsurf = linspace(min(x),max(x));
% zsurf = interp1(x,z,-.01,"linear");

%plot result
figure(1)
plot(x,z)
title('Simulated Reachable Region, X-Z Plane')
xlabel('x')
ylabel('z')
grid on
axis([-.07 .07 0 .08])