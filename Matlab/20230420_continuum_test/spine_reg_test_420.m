% restart
close all; clear;

% spine reg points
% c1 = [-106.07 0 -106.07]';
% c2 = [106.07 0 106.07]';
% c3 = [-106.07 0 106.07]';
% c4 = [106.07 0 -106.07]';
% b1 = [-6 5 0]';
% b2 = [6 5 0]';
% b3 = [0 5 6]';
% b4 = [0 5 -6]';

c1 = [-106.07 106.07 0]';
c2 = [106.07 -106.07 0]';
c3 = [-106.07 -106.07 0]';
c4 = [106.07 106.07 0]';
b1 = [-6 0 5]';
b2 = [6 0 5]';
b3 = [0 -6 5]';
b4 = [0 6 5]';

regpt_spine = [c1 c2 c3 c4 b1 b2 b3 b4];

% load tip file
penprobe_tip = load('penprobe.tip')';

% load registration points from aurora and apply tip compensation
reg_tab = readtable('20230420_reg_002_pen.csv','Delimiter',',');
regpt_aurora = nan(size(regpt_spine));
for pt_idx = 1:size(reg_tab,1)
    q = table2array(reg_tab(pt_idx,4:7))';
    t = table2array(reg_tab(pt_idx,8:10))';
    regpt_aurora(:,pt_idx) = t + quatrotate(q',penprobe_tip')';
end

% register aurora to model
[~,TF_aurora_to_model,rmse] = rigid_align_svd(regpt_aurora,regpt_spine)

% load spine model
spine_stl = stlread('dual_helix_5rotations.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);  % not sure why we can pass an array of more than 3 cols for vertices here, but it is OK with homogeneous coords...


% load path data
path_tab = readtable('20230420_spinecap_002.csv','Delimiter',',');
path_pts = [path_tab.Var3 path_tab.Var8 path_tab.Var9 path_tab.Var10]';
c = 1;
for i = 1:size(path_pts,2)
    if path_pts(1,i) == 1
        path_pts_tip(:,c) = path_pts(2:end,i);
        c = c+1;
    end
end
path_pts_model = hTF(path_pts_tip,TF_aurora_to_model,0);

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


