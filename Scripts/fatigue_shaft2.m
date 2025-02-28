close all; clear; clc;

% TODO: Konverter inpit S_ut fra MPa til ksi og velg nærmeste verdi fra
% tabell, konverter så tilbake til MPa

% Input parameters
n_f = 2; % Safety factor
S_ut = 1080; % [MPa] Material ultimate tensilte strength PLACEHLDER VALUE
d_shaft = 130.8; % [mm] Shaft diameter PLACEHLDER VALUE
r_fillet = 1; % [mm] Fillet radius PLACEHLDER VALUE
D_d = 1.2; % PLACEHLDER VALUE
load_type = "Complex axial";  % ("Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial");
surface_finish = "Machined"; % ("Ground" "Machined" "Hot-rolled" "As-forged") For other types, see Machine Design page 368, Figure 6-26
reliability = 99; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)
operating_temperature = 22; % Celsius


% From mecOfMaterials_shaft2.m 
% (Values are not accurate and must be updated)
M_y_max = 1200*1e3; % [Nmm]
M_y_min = -1200*1e3; % [Nmm]
M_z_max = 3800*1e3; % [Nmm]
M_z_min = -3800*1e3; % [Nmm]

T_max = 22700*1e3; % [Nmm]
T_min =  22700*1e3; % [Nmm]

% Lecture 3 slide 7
M_max = sqrt(M_y_max^2 + M_z_max^2); % [Nmm]
M_min = -sqrt(M_y_min^2 + M_z_min^2); % [Nmm]
M_mean = (M_max + M_min)/2; % [Nmm]
M_amp = (M_max - M_min)/2; % [Nmm]

T_mean = (T_max + T_min)/2; % [Nmm]
T_amp = (T_max - T_min)/2; % [Nmm]




% Under er det mye feil, jobber med å rydde opp og tyde (lecture 5 slide 3)

% sigma_amp = K_f-bend*(32*M_amp)/(pi*d_shaft^3);
% sigma_mean = K_f-bend*(32*M_mean)/(pi*d_shaft^3);

% tau_amp = K_f-tor*(16*T_amp)/(pi*d_shaft^3);
% tau_mean = K_f-tor*(16*T_mean)/(pi*d_shaft^3);

% sigma_max = ?
% sigma_min = ?
% sigma_mean = (sigma_max + sigma_min)/2; % Mean (midrange) stess
% if sigma_mean ~= 0
%     fprintf("sigma_mean = ", sigma_mean, " Non nully reversed loading!")
% end
% sigma_amp = (sigma_max - sigma_min)/2; % Alternating stress amplitude


% Neubler's Constant for Steels (Machine Design, Table 6-6 page 382)
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240];% [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi = S_ut_ksi_table(1, 11); % [ksi], midlertidig "manuelt" utvalg, skal programeres til å finne nærmeste verdri fra MPa
a_sqrt_in = Neublers_table(S_ut_ksi); % [in^1/2]
a_sqrt_mm = a_sqrt_in * sqrt(25.4); % [mm^1/2] Usikker på om dette er rett. (Maskinelementer, lecture 3 slide 13) 

% Notch sensitivity factor (Maskinelementer, lecture 3 slide 13 & Machine Design, equation 6.13 page 381 )
q = 1/(1 + (a_sqrt_mm/sqrt(r_fillet))); 

%%% Stress concentration factors %%% (Maskinelementer, lecture 3 slide 12)
% Geometrical (theoretical) stress concentration factors:
% See appendix C for tables and values (side 1048-1049)
D_d_bend_key = [6.00, 3.00, 2.00, 1.50, 1.20, 1.10, 1.07, 1.05, 1.03, 1.02, 1.01];
D_d_bend_values = struct('A', num2cell([0.87868, 0.89334, 0.90879, 0.93836, 0.97098, 0.95120, 0.97527, 0.98137, 0.98061, 0.96048, 0.91938]),'b', num2cell([-0.33243, -0.30860, -0.28598, -0.25759, -0.21796, -0.23757, -0.20958, -0.19653, -0.18381, -0.17711, -0.17032]));
D_d_bend_table = dictionary(D_d_bend_key, D_d_bend_values);

A_bend = D_d_bend_table(D_d).A;
b_bend = D_d_bend_table(D_d).b;

D_d_tor_key = [2.00, 1.33, 1.20, 1.09];
D_d_tor_values = struct('A', num2cell([0.86331, 0.84897, 0.83425, 0.90337]),'b', num2cell([-0.23865, -0.23161, -0.21649, -0.12692]));
D_d_tor_table = dictionary(D_d_tor_key, D_d_tor_values);

A_tor = D_d_tor_table(D_d).A;
b_tor = D_d_tor_table(D_d).b;

D_d_axial_key = [2.00, 1.50, 1.30, 1.20, 1.15, 1.10, 1.07, 1.05, 1.02, 1.01];
D_d_axial_values = struct('A', num2cell([1.01470, 0.99957, 0.99682, 0.96272, 0.98084, 0.98450, 0.98498, 1.00480, 1.01220, 0.98413]),'b', num2cell([-0.30035, -0.28221, -0.25751, -0.25527, -0.22485, -0.20818, -0.19548, -0.17076, -0.12474, -0.10474]));
D_d_axial_table = dictionary(D_d_axial_key, D_d_axial_values);

A_axial = D_d_axial_table(D_d).A;
b_axial = D_d_axial_table(D_d).b;

K_t_bend = A_bend*(r_fillet/d_shaft)^b_bend;  % for normal stress, Appendix C
K_t_tor = A_tor*(r_fillet/d_shaft)^b_tor; % for shear stress
K_t_axial = A_axial*(r_fillet/d_shaft)^b_axial;

% Fatigue (dynamic) stress concentration foactors:
K_f = 1 + q * (K_t_bend - 1);      % for normal stress, Appendix C
K_fs = 1 + q * (K_t_tor - 1);    % for shear stress

%%% Correction factors %%%
% Load factor % (Maskinelementer, lecture 3 slide 43 & Machine Design, page 366)
C_load_table_key = ["Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial"];
C_load_table_value = [1 0.75 1 0.75 1];
C_load_table = dictionary(C_load_table_key, C_load_table_value);

C_load = C_load_table(load_type);

% Size factor (Maskinelementer, lecture 3 slide 44 & Machine Design, page 367)
if d_shaft <= 8
    C_size = 1;
elseif d_shaft <= 250
    C_size = 1.189*d_shaft^(-0.097);
else
    C_size = 0.6;
end

% Surface factor (Maskinelementer, lecture 3 slide 45 & Machine Design, page 369)
C_surf_table_key = ["Ground" "Machined" "Hot-rolled" "As-forged"];
C_surf_table_value = struct('A', num2cell([1.58 4.51 57.7 272]), 'b', num2cell([-0.085 -0.265 -0.718 -0.995]));
C_surf_table = dictionary(C_surf_table_key, C_surf_table_value);
C_surf_A = C_surf_table(surface_finish).A;
C_surf_b = C_surf_table(surface_finish).b;

C_surf = C_surf_A*S_ut^C_surf_b;
if C_surf > 1
    C_surf = 1;
end

% Temperature factor (Maskinelementer, lecture 3 slide 47 & Machine Design, page 371)
if operating_temperature <= 450
    C_temp = 1;
elseif operating_temperature <= 550
    C_temp = 1-0.0058*(operating_temperature-450);
else
    fprintf("Error: Operating temperatore too high\n")
end

% Reliability factor (Maskinelementer, lecture 3 slide 48 & Machine Design, page 371)
C_reliab_table_key = [50 90 95 99 99.9 99.99 99.999 99.9999];
C_reliab_table_value = [1 0.897 0.868 0.814 0.753 0.702 0.659 0.620];
C_reliab_table = dictionary(C_reliab_table_key, C_reliab_table_value);
C_reliab = C_reliab_table(reliability);

% Other factors
% See lecture 3 slide 49

% For steels with "knee" (Maskinelementer, lecture 3 slide 39)
if S_ut < 1400
    S_e_prime = 0.5*S_ut;
else
    S_e_prime = 700;
end

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % (Maskinelementer, lecture 4 slide 5)


d = ((16*n_f/pi)*(sqrt(4*(K_f*M_amp)^2+3*(K_fs*T_amp)^2)/S_e)+(sqrt(4*(K_f*M_mean)^2+3*(K_f*T_mean)^2/S_ut)))^(1/3)
% Stemmer det at K_f_bend = K_f og at K_f_tor = K_f_s ??

%% First itteration

% Estimating stress geometric concentration factors for preliminary stage
% Shoulder fillet sharp (r/d = 0.02, D/d = 1.5), lecture 5 slide 10, maskinelementer)
% K_t_bend = 2.7;
% K_t_tor = 2.2;
% K_t_axial = 3.0;

% Shoulder fillet well-rounded (r/d = 0.1, D/d = 1.5), lecture 5 slide 10, maskinelementer)
K_t_bend = 1.7;
K_t_tor = 1.5;
K_t_bend = 1.9;

% Fatigue concentration foactor, lecture 5 slide 11
% K_f_bend = 1 + q * (K_t_bend - 1); % q unknown
% K_f_tor = 1 + q * (K_t_tor - 1); % q unknown
% K_f_axial = 1 + q * (K_t_axial - 1); % q unknown

% Conservative estimate for preliminary stage
K_f_bend = K_t_bend;
K_f_tor = K_t_tor;
K_f_axial = K_t_bend;

% % Correction factors for preliminary stage, lecture 5 slide 12
% % S_e = C_load * C_size * C_surf * C_temp * C_reliab * S_ePrime;
C_load = 1; % Bending
C_size = 1;
% C_surf = C_surf;
% C_temp = C_temp; 
% C_reliab = C_reliab;

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % (Maskinelementer, lecture 4 slide 5)


% Lecture 5 slide 6
% sigma_rev = sigma_amp_prime/(1-(simga_mean_prime/S_ut)); % 
% n_f = S_e / sigma_rev; %

% First itteration
d1 = ((16*n_f/pi)*(sqrt(4*(K_f_bend*M_amp)^2+3*(K_f_tor*T_amp)^2)/S_e)+(sqrt(4*(K_f_bend*M_mean)^2+3*(K_f_tor*T_mean)^2/S_ut)))^(1/3)