clc;close all;clear;

%TODO:
% find actual datasheet

lifetime = 10; % [years]
work_cycle = 10; % [hours/day]

% teeth #
z_1 = 18;
z_2 = 79;
z_3 = 18;
z_4 = 71;

% gear ratios of stages
i_s1 = z_2/z_1;
i_s2 = z_4/z_3;
i_tot = i_s1 * i_s2;

% speed of gears [rpm]
n_1 = 1450; % input
n_2 = n_1 / i_s1;
n_3 = n_2;
n_4 = n_3/ i_s2;

% number of cycles through lifetime:
cycles_lifetime_1 = lifetime * 365.25 * work_cycle*60 * n_1
cycles_lifetime_2 = lifetime * 365.25 * work_cycle*60 * n_2
cycles_lifetime_4 = lifetime * 365.25 * work_cycle*60 * n_4

% Loads from mechOfMaterials_shaftx.m scripts:
% - radial [N]:
F_r_sh1 = 20;
F_r_sh2 = 20;
F_r_sh3 = 20;
% - axial  [N]:
F_a_sh1 = 5;
F_a_sh2 = 5;
F_a_sh3 = 5;

% Equivalent load P
% factors from figure 11-24 pg 705 machine design book, 20deg taper
V = 1.0; % rotation factor, rotating inner ring
X = 0.43; % radial factor
Y = 1; % thrust (axial) factor
e_f = 0.3; % axial/thrust factor minimum

P_sh1 = X * V * F_r_sh1 + Y * F_a_sh1 %eq 11.22a pg 704 machine design

% L_10 fatigue life
% for now using bearing size 6320 from figure 11-23 pg 702 machine design book
C_sh1 = 28500 * 4.4482216153; % lb -> N dynamic load rating
C_0_sh1 = 27000 * 4.4482216153; % lb -> N static load rating
K_R = 1; % Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design

L_10_sh1 = K_R*(C_sh1/P_sh1)^(10/3) % eq 11.20b, pg 701 machine design book
if cycles_lifetime_1 > L_10_sh1
    warning("bearings for shaft 1 has low lifetime")
end