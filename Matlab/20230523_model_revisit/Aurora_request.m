% restart
close all; clear; clc;

PROBE_ID_PENCOIL  = 0x3CC3F000;   % pen probe coil
PROBE_ID_SPINE    = 0x3D4C7400;   % Brook's spine coil

%read in aurora information
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);
spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);

% load spine model for plotting
spine_stl = stlread('dual_helix.STL');
[f_resamp,v_resamp] = reducepatch(spine_stl.ConnectivityList,spine_stl.Points,0.1);

fser = serialport('/dev/tty.usbserial-1120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);

pkt = [];
while(isempty(pkt))
    requestAuroraPacket(fser,[PROBE_ID_PENCOIL PROBE_ID_SPINE]);
    pkt = getAuroraPacket(fser,0.2);
end

clear fser;

assert(mod(length(pkt)-1,9) == 0,'Incorrect packet size!');
num_tforms = (length(pkt)-1)/9;
all_tfs = [];
for tf_idx = 1:num_tforms
    all_tfs(tf_idx).sn = pkt(2+9*(tf_idx-1));
    all_tfs(tf_idx).T_coil_to_aurora = eye(4);
    all_tfs(tf_idx).T_coil_to_aurora(1:3,1:3) = quat2matrix(pkt((3:6)+9*(tf_idx-1)));
    all_tfs(tf_idx).T_coil_to_aurora(1:3,4) = pkt((7:9)+9*(tf_idx-1));
    all_tfs(tf_idx).error = pkt(10+9*(tf_idx-1));
end
if(num_tforms)
    all_tfs.T_coil_to_aurora
else
    warning('No transforms returned!');
end

spinedata_in_auroraspace = zeros(3,1);

for i = 1:length(all_tfs)
    if all_tfs(i).sn == PROBE_ID_SPINE
        spinedata_in_auroraspace = (all_tfs(i).T_coil_to_aurora(1:3,4) + all_tfs(i).T_coil_to_aurora(1:3,1:3)*spinetip_in_coilspace);
    end
end

spinedata_in_modelspace = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);

figure(1)
plot3(spinedata_in_modelspace(1),spinedata_in_modelspace(2),spinedata_in_modelspace(3),'g.','MarkerSize',50)
hold on;
patch('Vertices',v_resamp,'Faces',f_resamp,'EdgeColor','k','FaceColor',"#0072BD",'LineWidth',0.01);