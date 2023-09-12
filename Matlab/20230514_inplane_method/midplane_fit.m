clear; close all;

%physical characteristics
L = 56; % (mm) from spine architecture 
d = 4; % (mm) "

%model Rc range
res =100;
minR = log10(2.5*L/(pi));
maxR = log10(1000);
Rc = -[logspace(minR,maxR,res)]';
phi = L./Rc;

%model x and z
z = Rc.*sin(L./Rc);
x = -Rc.*(1-cos(L./Rc));
dl = phi.*(Rc-d)-L;

%model fits
[zmod_dl, gofz] = fit(z,Rc,'poly3')
[xmod_dl, gofx] = fit(x,Rc,'poly4')



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
plot(zmod_dl,z,dl)
hold on; grid on;
%plot(zfit1_Rc,side1_pos(3,:),Rc_oneside)
%legend('model','model fit','data','data fit','Location','southeast')
ylabel('dl')
xlabel('z')

figure(2)
plot(xmod_dl,x,dl)
grid on;
ylabel('dl')
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