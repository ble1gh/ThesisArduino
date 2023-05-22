clear; close all;

%physical length
L = 56;

%model Rc range
res =100;
minR = log10(2*L/(pi));
maxR = log10(300);
Rc = -[logspace(log10(3000),maxR,res/2) logspace(maxR,minR,res)]';

%model x and z
z = Rc.*sin(L./Rc);
x = -Rc.*(1-cos(L./Rc));

%model fits
[zmod_Rc, gofz] = fit(z,Rc,'fourier4');
[xmod_Rc, gofx] = fit(x,Rc,'exp2');

%read in real data
xz_tab = readtable('transformed_output_0514_xz.csv','Delimiter',',');
xz_pos = table2array(xz_tab);
Rc_tab = readtable('Sweep_Rc.csv','Delimiter',',');
Rc_data = table2array(Rc_tab);

%separate the two sides of the sweep
side1_pos = -xz_pos(:,1:length(xz_pos)/2);
side1_pos(3,:) = -side1_pos(3,:);
side2_pos = xz_pos(:,length(xz_pos)/2+1:end);
Rc_oneside = Rc_data(length(Rc_data)/2+1:end);

%fit data
%[zfit1_Rc, gofz] = fit(side1_pos(3,:)',Rc_oneside','fourier5')
% [xfit1_Rc, gofx] = fit(side1_pos(1,10:end)',Rc_oneside(10:end)','exp2')
% [xfit2_Rc, gofx] = fit(side2_pos(1,10:end)',Rc_oneside(10:end)','exp2')

[xfit1_Rc, gofx] = fit(side1_pos(1,:)',Rc_oneside','exp2')
[xfit2_Rc, gofx] = fit(side2_pos(1,:)',Rc_oneside','exp2')

figure(1)
plot(zmod_Rc,z,Rc)
hold on; grid on;
%plot(zfit1_Rc,side1_pos(3,:),Rc_oneside)
%legend('model','model fit','data','data fit','Location','southeast')
ylabel('Rc')
xlabel('z')

figure(2)
plot(xfit1_Rc,side1_pos(1,:),Rc_oneside)
hold on; grid on;
plot(x,Rc,'Color',"#77AC30")
ylabel('Rc (mm)')
xlabel('x (mm)')
legend('data','data fit','model')
title('x vs Rc, Model vs Data')

figure(3)
plot(xfit1_Rc,side1_pos(1,:),Rc_oneside)
hold on; grid on;
plot(xfit2_Rc,side2_pos(1,:),Rc_oneside,'.')
ylabel('Rc')
xlabel('x')
legend('data fit 1','side 1','data fit 2','side 2','Location','east')