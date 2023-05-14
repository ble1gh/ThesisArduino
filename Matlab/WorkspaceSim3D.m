%Brook Leigh
%Simulation of 3D reachable area for one segment
clear

L = 56; % (m) from spine architecture 

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(3*L/(2*pi));
maxR = log10(200);
Rc = [logspace(minR,maxR,50)];

%2 alternative ways to get acceptable Rc values
% Rc = [linspace(3*L/(2*pi),2*L/pi,50) linspace(2*L/pi,.8,100) linspace(-.8,-2*L/pi,100) linspace(-2*L/pi,-3*L/(2*pi),50)]; 
% Rc = [linspace(3*L/(2*pi),2*L/pi,50) linspace(2*L/pi,.8,100)]; 

theta = linspace(0,2*pi,50);
phi = L./Rc;

%initiate position variables
x = zeros(length(theta),length(Rc));
y = zeros(length(theta),length(Rc));
z = zeros(length(theta),length(Rc));

%for a given theta, an arc is defined along the range of possible bend
%radiuses which each give an (x,y,z) coordinate
for i = 1:length(theta)
    x(i,:) = cos(theta(i))*Rc.*(1-cos(L./Rc));
    y(i,:) = sin(theta(i))*Rc.*(1-cos(L./Rc));
    z(i,:) = Rc.*sin(L./Rc);
end


%plot result
figure(1)
%surf(x,y,z,'FaceColor','interp','FaceAlpha',0.9,'EdgeColor','none')
surf(x,y,z,'FaceColor',"#4DBEEE",'FaceAlpha',0.7)
title('Simulated Reachable Region, 3D')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-70 70 -70 70 0 80])


