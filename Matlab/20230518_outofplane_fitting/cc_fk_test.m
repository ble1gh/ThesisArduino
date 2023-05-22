% constant curvature spine forward kinematics test

% restart
close all; clear; clc;

L = 65;  % [mm] overall spine length
d = 5;   % [mm] radius from center of spine to cable
dLa = 5; % [mm] displacement of cable that is positioned on +x axis, positive is pulling
dLb = 0; % [mm] displacement of cable that is positioned on +y axis, positive is pulling

% copmute phi, theta, r0
phi = (1/d)*sqrt(dLa^2 + dLb^2);
theta = atan2(dLb,dLa);
r0 = L/phi;

% compute how far other cables need to be let out
dLc = L-(r0 + d*cos(theta));
dLd = L-(r0 + d*sin(theta));

% compute tip transform
t = [r0*(1-cos(phi))*cos(theta); r0*(1-cos(phi))*sin(theta); r0*sin(phi)]   % [mm]
R1 = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
R2 = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
T = [R2*R1*R2' t; zeros(1,3) 1];

% compute spine arc for display (line along center of spine)
arc_ang = 0:0.01:phi;
arc_xyz = [r0*(1-cos(arc_ang)); zeros(size(arc_ang)); r0*(sin(arc_ang))];
arc_xyz = R2*arc_xyz;

% plot results
figure;
hold on; grid on; axis equal;
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
plotTriad(eye(4),10);
plotTriad(T,10);
plot3(arc_xyz(1,:),arc_xyz(2,:),arc_xyz(3,:),'-','LineWidth',3,'Color',[1 0 1]);
view([50 30]);