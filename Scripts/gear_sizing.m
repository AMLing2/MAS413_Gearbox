clc;close all;clear;

alpha = 20; % [deg] helix angle, psi
beta = 15;  % [deg] pressure angle, theta

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

c_dist_1 = 50e-3; % [m] center distance from first to second gear , temp val
c_dist_2 = 150e-3; % [m]


% m_1 = (d_g1 / z_1) * 1e3 % [mm] module for stage 1, calc with diameter d_g1
m_s1 = (2*c_dist_1 / (z_1 + z_2)) * 1e3; % [mm] module for stage 1
m_s2 = (2*c_dist_2 / (z_3 + z_4)) * 1e3; % [mm] module for stage 2

% pitch circles [mm]
d_g1 = m_s1 * z_1;
d_g2 = m_s1 * z_2;
d_g3 = m_s2 * z_3;
d_g4 = m_s2 * z_4;

% top (ht) and bottom (hf) heights [mm]
ht_1 = m_s1; 
ht_2 = m_s1; 
ht_3 = m_s2; 
ht_4 = m_s2; 
hf_1 = 1.25 * m_s1; 
hf_2 = 1.25 * m_s1; 
hf_3 = 1.25 * m_s2; 
hf_4 = 1.25 * m_s2; 

% addedum [mm]
dt_g1 = d_g1 + 2 * ht_1;
dt_g2 = d_g2 + 2 * ht_2;
dt_g3 = d_g3 + 2 * ht_3;
dt_g4 = d_g4 + 2 * ht_4;

% dedendum [mm]
df_g1 = d_g1 - 2 * hf_1;
df_g2 = d_g2 - 2 * hf_2;
df_g3 = d_g3 - 2 * hf_3;
df_g4 = d_g4 - 2 * hf_4;

%gearbox total length of gears
l_tot = (dt_g1 + d_g2/2 + d_g3/2 + dt_g4)/1e3
l_in_to_out = c_dist_1 + c_dist_2

% pitch [mm]
p_s1 = pi*m_s1;
p_s2 = pi*m_s2;

% tooth thickness [mm]
sn_s1 = p_s1/2 - 0.05 * m_s1;
sn_s2 = p_s2/2 - 0.05 * m_s2;

% hatch width [mm]
en_s1 = p_s1/2 + 0.05 * m_s1;
en_s2 = p_s2/2 + 0.05 * m_s2;

% Torques in each gear [Nmm]
T_1 = 82.32152229 * 1e3; % [Nm] -> [Nmm]
T_2 = T_1 * i_s1;
T_3 = T_2;
T_4 = T_3 * i_s2;

%% helical gear calcs

%helical modules
mt_s1 = m_s1 / cosd(alpha);
mt_s2 = m_s2 / cosd(alpha);

%width of helical gears
lambda = 10; % width factor, 8-12, pg 17 lec 1
b_s1 = mt_s1 * lambda;
b_s2 = mt_s2 * lambda;

% transverse pitch [mm] 
pt_s1 = p_s1 / cosd(alpha);
pt_s2 = p_s2 / cosd(alpha);

% axial pitch [mm] 
px_s1 = p_s1 / sind(alpha);
px_s2 = p_s2 / sind(alpha);

% diameteral pitch [mm]
dp_s1 = pi/pt_s1;
dp_s2 = pi/pt_s2;

%% Bending stress sigma_b from lectures
A = 5; % [m/s] operating factor, tab 2 pg 6 lec 4
K_a = 1.25; % external dynamic factor, electric motor with light shock, tab 1 pg 5 lec 4
V_b = 1.7; % safety factor, lec 4 pg 8

%gear 1
gamma_1 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
V_t_1 = n_1 * ((2*pi)/60) * ((d_g1*1e-3)/2); %pitch speed [m/s] 
K_V_1 = (A + V_t_1) / A; % dynamic factor
F_th_1 = T_1 / (d_g1 / 2); % theoretical tangential force component [N]
%gear 2
gamma_2 = 2.24; % teeth form factor, 79 teeth, tab 3 pg 7 lec 4
V_t_2 = n_2 * ((2*pi)/60) * ((d_g2*1e-3)/2); %pitch speed [m/s] 
K_V_2 = (A + V_t_2) / A; % dynamic factor
F_th_2 = T_2 / (d_g2 / 2); % theoretical tangential force component [N]
%gear 3
gamma_3 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
V_t_3 = n_3 * ((2*pi)/60) * ((d_g3*1e-3)/2); %pitch speed [m/s] 
K_V_3 = (A + V_t_3) / A; % dynamic factor
F_th_3 = T_3 / (d_g3 / 2); % theoretical tangential force component [N]
%gear 4
gamma_4 = (2.30 + 2.24) / 2; % teeth form factor, 71 teeth, tab 3 pg 7 lec 4
V_t_4 = n_4 * ((2*pi)/60) * ((d_g4*1e-3)/2); %pitch speed [m/s] 
K_V_4 = (A + V_t_4) / A; % dynamic factor
F_th_4 = T_4 / (d_g4 / 2); % theoretical tangential force component [N]

% Material bending limits
sigma_b_lim_mat_list = [160,210,220,250,300,310,410,410]; % [MPa]
sigma_b_lim_mat_names = ["Fe 430", "Fe 590", "C 45 N", "C 60 N",...
    "34 Cr 4 V", "42 CrMo 4 V", "16 MnCr 5", "15 CrNi 6"];
sigma_b_lim_mat = dictionary(sigma_b_lim_mat_names,sigma_b_lim_mat_list);
sigma_b_lim = sigma_b_lim_mat("Fe 590") / V_b;

% calculation of modules based on sigma_b_lim

mt_1 = sqrt((F_th_1 * K_a * K_V_1 * gamma_1)/(sigma_b_lim * lambda));
mt_2 = sqrt((F_th_2 * K_a * K_V_2 * gamma_2)/(sigma_b_lim * lambda));
mt_3 = sqrt((F_th_3 * K_a * K_V_3 * gamma_3)/(sigma_b_lim * lambda));
mt_4 = sqrt((F_th_4 * K_a * K_V_4 * gamma_4)/(sigma_b_lim * lambda));
mt_s1 = max([mt_1,mt_2])
mt_s2 = max([mt_3,mt_4])
return
%calculating bending stress with known module
sigma_b_1 = (F_th_1* K_a * K_V_1 * gamma_1)/(b_s1*mt_s1) % [MPa]


%% Bending stress sigma_b - from mechanics book
return
%TODO:
%pg 753 machine design - assumptions
%contact ratio between 1 and 2
% stages must obey tables 12-4 12-5

% gemoetry factor J (helical) pg 798 mechanical design
% Tab 13-2 (psi = 20, theta = 20)
J_1 = 0.47; % p = 17 teeth, g = 55 teeth
J_2 = 0.54;
J_3 = 0.47; % p = 17 teeth, g = 55 teeth
J_4 = 0.54;

% dynamic factor K_V
Q_v = 9; % 6 to 11, based on lowest quality gear in stage, 12.17b pg 754
K_V_1 = dynamicFactor(n_1,Q_v,d_g1);
K_V_2 = dynamicFactor(n_2,Q_v,d_g2);
K_V_3 = dynamicFactor(n_3,Q_v,d_g3);
K_V_4 = dynamicFactor(n_4,Q_v,d_g4);

%TODO: continue with factors in page 758 machine design +

return
% bending stress of gear, (12.15si) pg 753 (spur) and pg 796 (helical)
sigma_b_s1 = (W_t/(F*m*J)) * (k_a*K_m/K_V) * K_s * K_B * K_I;


%% FUNCTIONS
function [K_V] = dynamicFactor(n,Q_V,d_g)
% Q_V = 6 to 11, based on lowest quality gear in stage, 12.17b pg 754
    V_t = n * ((2*pi)/60) * ((d_g*1e-3)/2); %pitch speed m/s 
    B = ((12 - Q_V)^(2/3))/4;
    A = 50 + 56 * (1-B);
    K_V = (A/(A + sqrt(200 * V_t)))^B;
end