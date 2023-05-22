

%read in coefficient data
x_coeffs_Q1_tab = readtable("x_coeffs_Q1.csv");
x_coeffs_Q1 = table2array(x_coeffs_Q1_tab);
x_coeffs_Q2_tab = readtable("x_coeffs_Q2.csv");
x_coeffs_Q2 = table2array(x_coeffs_Q2_tab);
x_coeffs_Q3_tab = readtable("x_coeffs_Q3.csv");
x_coeffs_Q3 = table2array(x_coeffs_Q3_tab);
x_coeffs_Q4_tab = readtable("x_coeffs_Q4.csv");
x_coeffs_Q4 = table2array(x_coeffs_Q4_tab);

z_coeffs_Q1_tab = readtable("z_coeffs_Q1.csv");
z_coeffs_Q1 = table2array(z_coeffs_Q1_tab);
z_coeffs_Q2_tab = readtable("z_coeffs_Q2.csv");
z_coeffs_Q2 = table2array(z_coeffs_Q2_tab);
z_coeffs_Q3_tab = readtable("z_coeffs_Q3.csv");
z_coeffs_Q3 = table2array(z_coeffs_Q3_tab);
z_coeffs_Q4_tab = readtable("z_coeffs_Q4.csv");
z_coeffs_Q4 = table2array(z_coeffs_Q4_tab);

%aggregate to one array for x and z
x_coeffs = [x_coeffs_Q1(:,1:9) x_coeffs_Q2(:,1:9) x_coeffs_Q3(:,1:9) x_coeffs_Q4(:,1:10)];
z_coeffs = [z_coeffs_Q1(:,1:9) z_coeffs_Q2(:,1:9) z_coeffs_Q3(:,1:9) z_coeffs_Q4(:,1:10)];

%write coefficients to file
writematrix(x_coeffs,"x_coeffs_full.csv")
writematrix(z_coeffs,"z_coeffs_full.csv")


% points_record_Q1_tab = readtable("points_record_Q1.csv");
% points_record_Q1 = table2array(points_record_Q1_tab);
% points_record_Q2_tab = readtable("points_record_Q2.csv");
% points_record_Q2 = table2array(points_record_Q2_tab);
% points_record_Q3_tab = readtable("points_record_Q3.csv");
% points_record_Q3 = table2array(points_record_Q3_tab);
% points_record_Q4_tab = readtable("points_record_Q4.csv");
% points_record_Q4 = table2array(points_record_Q4_tab);
% 
% points_record_circle = [points_record_Q1 points_record_Q2 points_record_Q3 points_record_Q4];
% writematrix(points_record_circle,"points_record_circle.csv")

