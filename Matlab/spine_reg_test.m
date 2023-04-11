% restart
close all; clear; clc;

% spine reg points
regpt_spine = [-13.5 0 0; 0 -13.5 0; 13.5 0 0; 0 13.5 0]';

% load tip file
penprobe_tip = load('penprobe.tip')';

% load registration points from aurora and apply tip compensation
reg_tab = readtable('reg.csv','Delimiter',',');
regpt_aurora = nan(size(regpt_spine));
for pt_idx = 1:size(reg_tab,1)
    q = table2array(reg_tab(pt_idx,4:7))';
    t = table2array(reg_tab(pt_idx,8:10))';
    regpt_aurora(:,pt_idx) = t + quatrotate(q,penprobe_tip);
end

% register aurora to model
[~,TF_aurora_to_model,rmse] = rigid_align_svd(regpt_aurora,regpt_spine)

% load spine model
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);  % not sure why we can pass an array of more than 3 cols for vertices here, but it is OK with homogeneous coords...


% load path data
path_tab = readtable('path.csv','Delimiter',',');
path_pts = [path_tab.Var8 path_tab.Var9 path_tab.Var10]';
path_pts_model = hTF(path_pts,TF_aurora_to_model,0);

% plot
figure;
hold on; grid on; axis equal;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor','g','LineWidth',0.01);
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
plot3(path_pts_model(1,:),path_pts_model(2,:),path_pts_model(3,:),'.')
view([53,32]);

% figure;
% hold on; grid on; axis equal;
% plot3(path_pts_model(1,:),path_pts_model(2,:),path_pts_model(1,:),'.')


