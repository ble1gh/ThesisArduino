%Brook Leigh
%Simulation of 3D reachable area for one segment
clear; close all;

L = 56; % (m) from spine architecture 

%These values give a turn angle (theta) and a bend amount (Rc)
minR = log10(2.5*L/(pi));
maxR = log10(1000);
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
title('Simulated Reachable Region, 3D Constant Curvature')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%axis([-70 70 -70 70 0 80])
hold on; grid on; axis equal;

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);
%spine tip in model space
spine_tip_modelspace = [0,0,56]';

figure(1)
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);

