clc;clear;close all;
n_1 = 1450; % Speed of input shaft 1 [rpm]
P_1 = 12.5e3; % Effect on shaft [W]
i_tot = 17.3; % total gear ratio
alpha = 20; % pressure angle [deg]
beta = 15; % helix angle [deg]

n_out = n_1/i_tot;


%calculations for internal gears
i_1 = 4.5; % first gear ratio of 2 stage from table 15-38, Lec2 pg12
i_2 = 17.3/i_1

z_1_min = 18; % minimum 18-20 from teacher
z_2 = i_1 * z_1_min;

if gcd(z_1_min,z_2) > 1
    gcd(z_1_min,z_2)
    warning("stage 1 gears are not relative prime")
end