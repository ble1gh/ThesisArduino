clear; close all;

%physical length
L = 56;

%model Rc range
res =100;
minR = log10(2*L/(pi));
maxR = log10(300);
Rc = -[logspace(minR,maxR,res/2)]';

%model x and z
z = Rc.*sin(L./Rc);
x = -Rc.*(1-cos(L./Rc));

%model fits
[zmod_Rc, gofz] = fit(z,Rc,'fourier4');
[xmod_Rc, gofx] = fit(x,Rc,'exp2');

%read in real data
xz_tab = readtable('transformed_output_0514_midplane.csv','Delimiter',',');
xz_pos = table2array(xz_tab);
Rc_tab = readtable('Sweep_Rc.csv','Delimiter',',');
Rc_data = table2array(Rc_tab);

%separate the two sides of the sweep
side1_pos = -xz_pos(:,1:length(xz_pos)/2);
side1_pos(3,:) = -side1_pos(3,:);
side2_pos = xz_pos(:,length(xz_pos)/2+1:end);
Rc_oneside = Rc_data(length(Rc_data)/2+1:end);

%create combined data set from both theoretically identical sweeps
both_sides = [side1_pos side2_pos];
Rc_bothside = [Rc_oneside Rc_oneside];

%fit data
%[zfit1_Rc, gofz] = fit(side1_pos(3,:)',Rc_oneside','fourier5')
[xfit1_Rc, gofx] = fit(side1_pos(1,:)',Rc_oneside','exp2')
[xfit2_Rc, gofx] = fit(side2_pos(1,:)',Rc_oneside','exp2')
[xfitboth_Rc, gofx] = fit(both_sides(1,:)',Rc_bothside','exp2')

figure(1)
plot(zmod_Rc,z,Rc)
hold on; grid on;
%plot(zfit1_Rc,side1_pos(3,:),Rc_oneside)
%legend('model','model fit','data','data fit','Location','southeast')
ylabel('Rc')
xlabel('z')

figure(2)
plot(xmod_Rc,x,Rc)
hold on; grid on;
plot(side1_pos(1,:),Rc_oneside,'.',side2_pos(1,:),Rc_oneside,'.')
ylabel('Rc')
xlabel('x')
legend('model','model fit','data side 1','data side 2')

figure(3)
plot(xfit1_Rc,side1_pos(1,:),Rc_oneside)
hold on; grid on;
plot(xfit2_Rc,side2_pos(1,:),Rc_oneside,'.')
ylabel('Rc')
xlabel('x')
legend('data fit 1','side 1','data fit 2','side 2','Location','east')

figure(4)
plot(xfitboth_Rc,both_sides(1,:),Rc_bothside)
hold on; grid on;
ylabel('Rc')
xlabel('x')
%legend('data fit 1','side 1','data fit 2','side 2','Location','east')