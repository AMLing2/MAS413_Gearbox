close all; clear; clc;

% Input parameters
S_ut = 2000; % [MPa] Placeholder value
d = 20; % [mm] Placeholder shaft diameter value
surface_finish = "Ground"; % ("Ground" "Machined" "Hot-rolled" "As-forged") For other types, see Machine Design page 368, Figure 6-26
operating_temperature = 160; % Celsius
reliability = 99; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)


% From mecOfMaterials_shaft2.m

% sigma_max = ?
% sigma_min = ?
% sigma_mean = (sigma_max + sigma_min)/2; % Mean (midrange) stess 
% sigma_amp = (sigma_max - sigma_min)/2; % Alternating stress amplitude
% 
% M_y = ?
% M_z = ?
% 
% M = sqrt(M_y^2 + M_z^2);

% Neubler's Constant for Steels (Machine Design, Table 6-6 page 382)
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240];% [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi = S_ut_ksi_table(1, 1); % [ksi]
a_sqrt_in = Neublers_table(S_ut_ksi); % [in^1/2]
a_sqrt_mm = a_sqrt_in * sqrt(25.4); % [mm^1/2] Usikker på om dette er rett. (Maskinelementer, lecture 3 slide 13) 

% Notch sensitivity factor % (Maskinelementer, lecture 3 slide 13 & Machine Design, equation 6.13 page 381 )
% q = 1/(1 + (a_sqrt_mm/sqrt(r))); 

%%% Stress concentration factors %%% (Maskinelementer, lecture 3 slide 12)
% Geometrical (theoretical) stress concentration factors:
% K_t = 1;  % for normal stress, Appendix C
% K_ts = 1; % for shear stress

% Fatigue (dynamic) stress concentration foactors:
% K_f = 1 + q * (K_t - 1);      % for normal stress, Appendix C
% K_fs = 1 + q * (K_ts - 1);    % for shear stress

%%% Correction factors %%%
% Load factor % (Maskinelementer, lecture 3 slide 43 & Machine Design, page 366)
C_load_table_key = ["Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial"];
C_load_table_value = [1 0.75 1 0.75 1];
C_load_table = dictionary(C_load_table_key, C_load_table_value);

% Size factor (Maskinelementer, lecture 3 slide 44 & Machine Design, page 367)
if d <= 8
    C_size = 1;
elseif d <= 250
    C_size = 1.189*d^(-0.097);
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

%%%

% Estimating stress geometric concentration factors for preliminary stage
% Shoulder fillet sharp (r/d = 0.02, D/d = 1.5), lecture 5 slide 10, maskinelementer)
% K_tBend = 2.7;
% K_tTor = 2.2;
% K_tAxial = 3.0;

% Fatigue concentration foactor, lecture 5 slide 11
% K_fBend = 1 + q * (K_tBend - 1); % q unknown
% K_fTor = 1 + q * (K_tTor - 1); % q unknown
% K_fAxial = 1 + q * (K_tAxial - 1); % q unknown
% 
% Conservative estimate for preliminary stage
% K_fBend = K_tBend;
% K_fTor = K_tTor;
% K_fAxial = K_tAxial;
% 
% % Correction factors for preliminary stage, lecture 5 slide 12
% % S_e = C_load * C_size * C_surf * C_temp * C_reliab * S_ePrime;
% C_load = 1; % Bending
% C_size = 1;
% % C_surf = A * S_ut^b; % A and b from table 6-3 page 369
% % C_temp = ?
% % C_reliab = ?
% 
