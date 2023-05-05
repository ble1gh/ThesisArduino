close all; clear all;

c1 = [-106.07 106.07 -5]';
c2 = [106.07 -106.07 -5]';
c3 = [-106.07 -106.07 -5]';
c4 = [106.07 106.07 -5]';
b1 = [-6 0 0]';
b2 = [6 0 0]';
b3 = [0 -6 0]';
b4 = [0 6 0]';

regpt_spine = [c1 c2 c3 c4 b1 b2 b3 b4];

%% 

% load tip file
penprobe_tip = load('tipcal.tip');

% load registration points from aurora and apply tip compensation
reg_tab = readtable('20230427_reg_001.csv','Delimiter',',');
regpt_aurora = nan(size(regpt_spine));
for pt_idx = 1:size(reg_tab,1)
    q = table2array(reg_tab(pt_idx,4:7));
    t = table2array(reg_tab(pt_idx,8:10));
    regpt_aurora(:,pt_idx) = t + quatrotate(q,penprobe_tip);
end

% register aurora to model
[~,TF_aurora_to_model,rmse] = rigid_align_svd(regpt_aurora(:,2:4:32),regpt_spine)

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);  % not sure why we can pass an array of more than 3 cols for vertices here, but it is OK with homogeneous coords...

%spine tip in model space
spine_tip_modelspace = [0,0,56]';

%coil to aurora transform, this is just using the first coil point, could
%filter to use more accurate
qcoil = table2array(reg_tab(1,4:7));
tcoil = table2array(reg_tab(1,8:10));
TF_coil_to_aurora = eye(4);
TF_coil_to_aurora(1:3,1:3) = quat2matrix(qcoil);
TF_coil_to_aurora(1:3,4) = tcoil;

%model to coil and find the spine tip in coil space
TF_model_to_coil = inv(TF_coil_to_aurora)*inv(TF_aurora_to_model);
spinetip_in_coilspace = hTF(spine_tip_modelspace,TF_model_to_coil,0)

coil_at_registration = tcoil + quatrotate(qcoil,spinetip_in_coilspace');

%for coil tracking data, apply spine tip offset
rawcoil_tab = readtable('output.csv','Delimiter',',');
spinedata_in_auroraspace = zeros(3,size(rawcoil_tab,1));
for coil_idx = 1:size(rawcoil_tab,1)
    q = table2array(rawcoil_tab(coil_idx,3:6));
    t = table2array(rawcoil_tab(coil_idx,7:9));
    spinedata_in_auroraspace(:,coil_idx) = (t + quatrotate(q,spinetip_in_coilspace'))';
end

%transform spine data from aurora space to model space
spinedata_in_modelspace = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);
tip_at_registration_modelspace = hTF(coil_at_registration',TF_aurora_to_model)
spinedata_in_modelspace = [spinedata_in_modelspace tip_at_registration_modelspace];

%plot
figure(1)
plot3(spinedata_in_modelspace(1,:),spinedata_in_modelspace(2,:),spinedata_in_modelspace(3,:),'.','Color',"#77AC30");
hold on; grid on; axis equal;
plot3(tip_at_registration_modelspace(1),tip_at_registration_modelspace(2),tip_at_registration_modelspace(3),'g.','MarkerSize',50);
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);
xlabel('x')
ylabel('y')
zlabel('z')

%write results to file to save
writematrix(spinedata_in_modelspace,'transformed_output_0427.csv')

hold off;
