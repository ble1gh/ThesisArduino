clear;

%physical system parameters
L = 56; % (mm) from spine architecture 
d = 4;  % (mm) "
disk_diameter = 15; % (mm) diameter of the disk
PWMrange = 500-100;
setmid = 350;

PROBE_ID_PENCOIL  = 0x3CC3F000;   % pen probe coil
PROBE_ID_SPINE    = 0x3D4C7400;   % Brook's spine coil

%read in coefficients and workspace
x_coeffs_tab = readtable("x_coeffs_full.csv");
z_coeffs_tab = readtable("zfits_poly3_full_circle.csv");
x_coeffs = table2array(x_coeffs_tab);
z_coeffs = table2array(z_coeffs_tab);
workspace_tab = readtable("known_surface.csv");
workspace = table2array(workspace_tab);

%read in aurora information
TF_aurora_to_model_tab = readtable('TF_aurora_to_model.csv');
TF_aurora_to_model = table2array(TF_aurora_to_model_tab);
spinetip_in_coilspace_tab = readtable('spinetip_in_coilspace.csv');
spinetip_in_coilspace = table2array(spinetip_in_coilspace_tab);
% load tip file
penprobe_tip = load('tipcal.tip');

%open serial ports
fser = serialport('/dev/tty.usbserial-1120',115200,'DataBits',8,'FlowControl','none','StopBits',1,'Timeout',0.001);
device = serialport("/dev/tty.usbmodem11301",115200);
flush(device);

rep = 1;

while 1
    %find pen location
    pkt = [];
    while(isempty(pkt))
        requestAuroraPacket(fser,[PROBE_ID_PENCOIL PROBE_ID_SPINE]);
        pkt = getAuroraPacket(fser,0.2);
    end
    
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
    pendata_in_auroraspace = zeros(3,1);
    
    for i = 1:length(all_tfs)
        if all_tfs(i).sn == PROBE_ID_SPINE
            spinedata_in_auroraspace = (all_tfs(i).T_coil_to_aurora(1:3,4) + all_tfs(i).T_coil_to_aurora(1:3,1:3)*spinetip_in_coilspace);
        else
            pendata_in_auroraspace = (all_tfs(i).T_coil_to_aurora(1:3,4) + all_tfs(i).T_coil_to_aurora(1:3,1:3)*penprobe_tip');
        end
    end
    
    spinedata_in_modelspace(:,rep) = hTF(spinedata_in_auroraspace,TF_aurora_to_model,0);
    pendata_in_modelspace(:,rep) = hTF(pendata_in_auroraspace,TF_aurora_to_model,0)
    
    
    %create point requests and conver to theta and r
    pen = pendata_in_modelspace(:,rep);
    dist = vecnorm(workspace-pen);
    [val,indx] = min(dist);
    req(:,rep) = workspace(:,indx)
    r_req = sqrt(req(1,rep)^2+req(2,rep)^2)

    if r_req > 33
        r_req = 33
    end

    theta = atan2(req(2,rep),req(1,rep))
    if theta < 0
        theta = 2*pi+theta
    end
    
    %linearly interpolate between known theta curves
    ax = interp1(x_coeffs(1,:),x_coeffs(2,:),theta,"linear");
    bx = interp1(x_coeffs(1,:),x_coeffs(3,:),theta,"linear");
    cx = interp1(x_coeffs(1,:),x_coeffs(4,:),theta,"linear");
    dx = interp1(x_coeffs(1,:),x_coeffs(5,:),theta,"linear");
    ex = interp1(x_coeffs(1,:),x_coeffs(6,:),theta,"linear");
    az = interp1(z_coeffs(1,:),z_coeffs(2,:),theta,"linear");
    bz = interp1(z_coeffs(1,:),z_coeffs(3,:),theta,"linear");
    cz = interp1(z_coeffs(1,:),z_coeffs(4,:),theta,"linear");
    dz = interp1(z_coeffs(1,:),z_coeffs(5,:),theta,"linear");
    
    %create cable change array
    dlplane(rep) = ax*r_req.^4 + bx*r_req.^3 + cx*r_req.^2 + dx*r_req.^1 + ex;
    dl(1:4,rep) = [cos(theta)*dlplane(rep); -cos(theta)*dlplane(rep); sin(theta)*dlplane(rep); -sin(theta)*dlplane(rep)];
    
    %convert to PWM motor inputs
    dtheta = dl(1:4,rep)./(disk_diameter/2);
    dPWM = -dtheta*PWMrange/pi;
    u(1:4,rep) = setmid;
    u(1:4,rep) = u(1:4,rep)+dPWM;

    write(device,u(:,rep),"uint16")
    count = size(u(:,rep),1);
    response = read(device,count,"uint16")

    % c1 = [50*cos(-pi/4) 50*sin(-pi/4) -5-1.57]';
    % c2 = [50*cos(pi/4) 50*sin(pi/4) -5-1.57]';
    % c3 = [50*cos(3*pi/4) 50*sin(3*pi/4) -5-1.57]';
    % c4 = [50*cos(5*pi/4)+10*cos(-pi/4) 50*sin(5*pi/4)+10*sin(-pi/4) -5-1.57]';
    % 
    % regpt_spine = [c1 c2 c3 c4];
    % 
    % vecnorm(regpt_spine-pendata_in_modelspace(:,rep))

    pause(0.1)
    
    rep = rep + 1;
end

clear fser; clear device;