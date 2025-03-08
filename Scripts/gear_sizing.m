clc;close all;clear;

%TODO: 
% contact stress, usually sigma_o > sigma_b - done
% change m_n to m_t for all calcs - done
% check i_tot is 1% of 17.3 - done
% contact ratio between 1 and 2 - checked
% material factor see line 73
% stages must obey tables 12-4 12-5 pg 737, add check z2,z4 < 1309 ?
% increase lambda to 14?

%%% Chosen parameters
material = "15 CrNi 6";
lambda = 12; % width factor, processed:  8-12, pg 17 lec 1

%from requirements:
beta = 15;   % [deg] helix angle, psi
alpha = 20;  % [deg] pressure angle, theta
P1 = 12.5e3; % [W] input power
i_tot_og = 17.3;
V_b = 1.7; % bending safety factor, lec 4 pg 8
V_o = 1.25; % contact safety factor, lec 4 pg 10

CR_min = 1.4; % machine design pg 738
a_CR_min = 1.15; % machine design pg 796

% teeth #
z_1 = 18;
z_2 = 79;
z_3 = 18;
z_4 = 71;

% gear ratios of stages
i_s1 = z_2/z_1;
i_s2 = z_4/z_3;
i_tot = i_s1 * i_s2;

% check if new i_tot is within 1% of requirement
if ( abs(i_tot - i_tot_og)/i_tot_og ) > 0.01
    error("i_tot is greater than 1% of requirement")
end

% speed of gears [rpm]
n_1 = 1450; % input
n_2 = n_1 / i_s1;
n_3 = n_2;
n_4 = n_3/ i_s2;

% Torques in each gear [Nmm]
T_1 = 82.32152229 * 1e3; % [Nm] -> [Nmm]
T_2 = T_1 * i_s1;
T_3 = T_2;
T_4 = T_3 * i_s2;

%%%% add helical module calcs

% module of elasticity and material standards (no price):
% Fe 430: https://matweb.com/search/DataSheet.aspx?MatGUID=c2ba59bb365942a7b6da46f1cee370b8
% Fe 590: https://matweb.com/search/DataSheet.aspx?MatGUID=1dc0414bd1ea4061a5dc09382c455e2a
% C 45 N: https://matweb.com/search/DataSheet.aspx?MatGUID=2ca9b42e83894e8a8a61385fd7da63ae
% C 60 N: https://matweb.com/search/DataSheet.aspx?MatGUID=0a471605c1324daa910855e54a21fab3
% 34 Cr 4 V: https://matweb.com/search/DataSheet.aspx?MatGUID=4877d405464f448a96786c8cbd00d3b5
% 42 CrMo 4 V: https://matweb.com/search/DataSheet.aspx?MatGUID=38108bfd64c44b4c9c6a02af78d5b6c6
% 16 MnCr 5: https://matweb.com/search/DataSheet.aspx?MatGUID=2ab813ffa05d40329dffe0ee7f58b5de
% 15 CrNi 6: https://matweb.com/search/DataSheet.aspx?MatGUID=9ab3bf332758468ab36010790bd94349 ?

% Material limits, table 5 forelesning 4, pg 11
sigma_b_lim_mat_list = [160,210,220,250,300,310,410,410]; % [MPa]
sigma_o_lim_mat_list = [430,520,540,610,715,760,1600,1900]; % [MPa]
E_mat_list = [200 200 210 210 205 205 200 210]; %GPa
mat_names = ["Fe 430", "Fe 590", "C 45 N", "C 60 N",...
    "34 Cr 4 V", "42 CrMo 4 V", "16 MnCr 5", "15 CrNi 6"];
sigma_b_lim_mat = dictionary(mat_names,sigma_b_lim_mat_list);
sigma_o_lim_mat = dictionary(mat_names,sigma_o_lim_mat_list);
E_mat_dic = dictionary(mat_names,E_mat_list);

K_L = 1.0; % lubrication factor pg 10 lec 4

sigma_b_lim = sigma_b_lim_mat(material) / V_b;
% Z_v is calculated in the module_calc function
sigma_o_lim = (sigma_o_lim_mat(material) / V_o) * K_L;

%% Bending stress sigma_b from lectures
A = 5; % [m/s] operating factor, tab 2 pg 6 lec 4
K_a = 1.25; % external dynamic factor, electric motor with light shock, tab 1 pg 5 lec 4
F_w = sqrt(0.35 * (E_mat_dic(material)*1e3)); % [sqrt(N/mm^2)] material factor, lec 4 pg 9
F_c = 1.76; % edge form factor for alpha = 20, lec 4 pg 9

% gammas for gears
gamma_1 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
gamma_2 = 2.24; % teeth form factor, 79 teeth, tab 3 pg 7 lec 4
gamma_3 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
gamma_4 = (2.30 + 2.24) / 2; % teeth form factor, 71 teeth, tab 3 pg 7 lec 4

% calculating normal modules for each gear
m_n_1 = module_calc(0.01, sigma_b_lim, sigma_o_lim, z_1, n_1, T_1, A, ...
        K_a, lambda, gamma_1, F_w, F_c, "pinion", i_s1);
m_n_2 = module_calc(0.01, sigma_b_lim, sigma_o_lim, z_2, n_2, T_2, A, ...
        K_a, lambda, gamma_2, F_w, F_c, m_n_1 * z_1, i_s1);
m_n_3 = module_calc(0.01, sigma_b_lim, sigma_o_lim, z_3, n_3, T_3, A, ...
        K_a, lambda, gamma_3, F_w, F_c, "pinion", i_s2);
m_n_4 = module_calc(0.01, sigma_b_lim, sigma_o_lim, z_4, n_4, T_4, A, ...
        K_a, lambda, gamma_4, F_w, F_c, m_n_3 * z_3, i_s2);

% converting to transverse (helical) module
mt_1 = m_n_1 / cosd(beta);
mt_2 = m_n_2 / cosd(beta);
mt_3 = m_n_3 / cosd(beta);
mt_4 = m_n_4 / cosd(beta);
mt_s1 = max([mt_1,mt_2]) % 3.1180 w/ 15 CrNi 6
mt_s2 = max([mt_3,mt_4]) % 4.5334 w/ 15 CrNi 6

%%%%%%%% sizing calcs for helical gears


% pitch circle diameters [mm]
d_g1 = mt_s1 * z_1;
d_g2 = mt_s1 * z_2;
d_g3 = mt_s2 * z_3;
d_g4 = mt_s2 * z_4;

% top (ht) and bottom (hf) heights [mm]
ht_1 = mt_s1;
ht_2 = mt_s1; % see figure 12-8 pg 734 machine element
ht_3 = mt_s2; 
ht_4 = mt_s2; 
hf_1 = 1.25 * mt_s1;
hf_2 = 1.25 * mt_s1; 
hf_3 = 1.25 * mt_s2; 
hf_4 = 1.25 * mt_s2; 

% addedum circle [mm]
dt_g1 = d_g1 + 2 * ht_1;
dt_g2 = d_g2 + 2 * ht_2;
dt_g3 = d_g3 + 2 * ht_3;
dt_g4 = d_g4 + 2 * ht_4;

% dedendum circle [mm]
df_g1 = d_g1 - 2 * hf_1;
df_g2 = d_g2 - 2 * hf_2;
df_g3 = d_g3 - 2 * hf_3;
df_g4 = d_g4 - 2 * hf_4;

%gearbox total length of gears
l_tot = (dt_g1 + d_g2/2 + d_g3/2 + dt_g4)/1e3; % [m]

% pitch [mm]
p_s1 = pi*mt_s1;
p_s2 = pi*mt_s2;

% tooth thickness [mm]
sn_s1 = p_s1/2 - 0.05 * mt_s1;
sn_s2 = p_s2/2 - 0.05 * mt_s2;

% hatch width [mm]
en_s1 = p_s1/2 + 0.05 * mt_s1;
en_s2 = p_s2/2 + 0.05 * mt_s2;

% axial pitch [mm] 
px_s1 = p_s1 / sind(beta);
px_s2 = p_s2 / sind(beta);

% diameteral pitch [mm]
dp_s1 = pi/p_s1;
dp_s2 = pi/p_s2;

% width of helical gears [mm]
b_s1 = mt_s1 * lambda;
b_s2 = mt_s2 * lambda;

% rough sum of material volume for gears [cm^3]
material_sum = (pi * b_s1 * (((d_g1/2)^2)+(d_g2/2)^2) + ...
               pi * b_s2 * (((d_g3/2)^2)+(d_g4/2)^2)) * 1e-3
material_price = 4.3175e3; % [dollar/m^3]
rough_mat_price = material_sum*1e-6 * material_price

% fillet radius
pd_in_s1 = 25.4/mt_s1; % [1/in] machine design eq 12.4d pg735
r_f_1 = (0.3/pd_in_s1)*25.4; %[in] -> [mm] machine design tab 12-1 pg735
pd_in_s2 = 25.4/mt_s2; % [1/in] machine design eq 12.4d pg735
r_f_2 = (0.3/pd_in_s2)*25.4; %[in] -> [mm] machine design tab 12-1 pg735

% center distance
C_s1 = d_g1/2 + d_g2/2; % [mm]
C_s2 = d_g3/2 + d_g4/2; % [mm]

%%% contact ratio %%%
Z_s1 =  sqrt((d_g1/2 + ht_1)^2 - (d_g1/2 * cosd(alpha))^2) + ...
        sqrt((d_g2/2 + ht_2)^2 - (d_g2/2 * cosd(alpha))^2) ...
        - C_s1 * sind(alpha); % eq 12.2 machine design pg 730
P_b_s1 = (pi * d_g1/ z_1) * cosd(alpha); % eq 12.3b machine design pg 734
CR_s1 = Z_s1/P_b_s1; % eq 12.7a machine design pg 738

Z_s2 =  sqrt((d_g3/2 + ht_3)^2 - (d_g3/2 * cosd(alpha))^2) + ...
        sqrt((d_g4/2 + ht_4)^2 - (d_g4/2 * cosd(alpha))^2) ...
        - C_s2 * sind(alpha); % eq 12.2 machine design pg 730
P_b_s1 = (pi * d_g3/ z_3) * cosd(alpha); % eq 12.3b machine design pg 734
CR_s2 = Z_s2/P_b_s1; % eq 12.7a machine design pg 738
if or( (CR_s1 < CR_min), (CR_s2 < CR_min) )
    error("Contact ratio is less than minimum")
end
%%% contact ratio %%%

%%% contact ratio - alternate from Lecture 3 page 12 %%%
    % More readable table from:
    % https://brainly.in/question/3177379
% Variable shift for consistency with table, copy existing
a_x = C_s1;                  % Center Distance
m_t = mt_s1;                 % Transverse Module
% alpha_t = atand( tand(alpha)/cosd(beta) );
alpha_t = alpha;             % Transverse (Helical) Angle
d_a1 = dt_g1;                % Addendum Diameter
d_a2 = dt_g2;                % Addendum Diameter
d_b1 = d_g1 * cosd(alpha_t); % Base Diameter
d_b2 = d_g2 * cosd(alpha_t); % Base Diameter
% From L3p12 - Common for all in table
alpha_wt = acosd( (d_b1 + d_b2) / (2*a_x) );
% From L3p12 - Contact Ratio of Helical Pair
epsilon_a = ( sqrt( (d_a1/2)^2 - (d_b1/2)^2 ) + ...
              sqrt( (d_a2/2)^2 - (d_b2/2)^2 ) - ...
              a_x * sind(alpha_wt) ) / ...
            pi*m_t*cosd(alpha_t);
%%% contact ratio - alternate from Lecture 3 page 12 %%%
contactRatioStep1 = table(CR_s1, epsilon_a)

%%% axial contact ratio %%%
m_f_s1 = (b_s1 * tand(beta)) / (pi * mt_s1); % eq 13.5 machine design pg 796
m_f_s2 = (b_s2 * tand(beta)) / (pi * mt_s2);
if or((m_f_s1 < 1),(m_f_s2 < 1))
    error("Axial contact ratio is under 1")
elseif or((m_f_s1 < a_CR_min),(m_f_s2 < a_CR_min))
    warning("Axial contact ratio is low")
end
%%% axial contact ratio %%%

% Display results
T = table([d_g1; d_g2; d_g3; d_g4], [b_s1; b_s1; b_s2; b_s2], ...
    'VariableNames', {'Diameter', 'Width'}, ...
    'RowNames', {'gear 1', 'gear 2', 'gear 3', 'gear 4'});

disp(T)

save("gear_sizes.mat")