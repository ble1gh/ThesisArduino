close all; clear all;

L = 55;

c1 = [50*cos(-pi/4) 50*sin(-pi/4) -5-1.57]';
c2 = [50*cos(pi/4) 50*sin(pi/4) -5-1.57]';
c3 = [50*cos(3*pi/4) 50*sin(3*pi/4) -5-1.57]';
c4 = [50*cos(5*pi/4)+10*cos(-pi/4) 50*sin(5*pi/4)+10*sin(-pi/4) -5-1.57]';

regpt_spine = [c1 c2 c3 c4];


% load tip file
penprobe_tip = load('tipcal.tip');

% load registration points from aurora and apply tip compensation
reg_tab = readtable('20230521_reg_001.csv','Delimiter',',');
regpt_aurora = nan(size(regpt_spine));
for pt_idx = 1:size(reg_tab,1)
    q = table2array(reg_tab(pt_idx,4:7));
    t = table2array(reg_tab(pt_idx,8:10));
    regpt_aurora(:,pt_idx) = t + quatrotate(q,penprobe_tip);
end

% register aurora to model
[~,TF_aurora_to_model,rmse] = rigid_align_svd(regpt_aurora(:,1:4:16),regpt_spine)

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);
%spine tip in model space
spine_tip_modelspace = [0,0,L]';

%coil to aurora transform, this is just using the first coil point, could
%filter to use more accurate
straight_tab = readtable('20230521_reg_001.csv','Delimiter',',');
qcoil = table2array(straight_tab(2,4:7));
tcoil = table2array(straight_tab(2,8:10));
TF_coil_to_aurora = eye(4);
TF_coil_to_aurora(1:3,1:3) = quat2matrix(qcoil);
TF_coil_to_aurora(1:3,4) = tcoil;

%model to coil and find the spine tip in coil space
TF_model_to_coil = inv(TF_coil_to_aurora)*inv(TF_aurora_to_model);
spinetip_in_coilspace = hTF(spine_tip_modelspace,TF_model_to_coil,0)

coil_at_registration = tcoil + quatrotate(qcoil,spinetip_in_coilspace');

writematrix(TF_aurora_to_model,'TF_aurora_to_model.csv')
writematrix(spinetip_in_coilspace,'spinetip_in_coilspace.csv')

%% 

%for coil tracking data, apply spine tip offset
rawcoil_tab = readtable('modelb_output.csv','Delimiter',',');
spinedata_in_auroraspace = zeros(3,size(rawcoil_tab,1));
for coil_idx = 1:size(rawcoil_tab,1)
    q = table2array(rawcoil_tab(coil_idx,3:6));
    t = table2array(rawcoil_tab(coil_idx,7:9));
    spinedata_in_auroraspace(:,coil_idx) = (t + quatrotate(q,spinetip_in_coilspace'))';
end

coil_data = table2array(rawcoil_tab);
coil_data_modelspace = hTF(coil_data(:,7:9)',TF_aurora_to_model);

%transform spine data from aurora space to model space
spinedata_in_modelspace = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);
tip_at_registration_modelspace = hTF(coil_at_registration',TF_aurora_to_model)
%spinedata_in_modelspace = [spinedata_in_modelspace tip_at_registration_modelspace];

regpt_aurora_modelspace = hTF(regpt_aurora,TF_aurora_to_model);

%plot
figure(1)
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','Color',"#77AC30");
hold on; grid on; axis equal;
%plot3(coil_data_modelspace(1,:),coil_data_modelspace(2,:),coil_data_modelspace(3,:))
%plot3(regpt_spine(1,:),regpt_spine(2,:),regpt_spine(3,:))
%plot3(regpt_aurora_modelspace(1,:),regpt_aurora_modelspace(2,:),regpt_aurora_modelspace(3,:),'.','MarkerSize',30)
%plot3(tip_at_registration_modelspace(1),tip_at_registration_modelspace(2),tip_at_registration_modelspace(3),'g.','MarkerSize',50);
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x')
ylabel('y')
zlabel('z')

%write results to file to save
writematrix(spinedata_in_modelspace,'transformed_output_0517_xz.csv')
%writematrix(TF_coil_to_aurora,'TF_coil_to_aurora.csv')


hold off;
