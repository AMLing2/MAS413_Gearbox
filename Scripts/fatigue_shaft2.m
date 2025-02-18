close all; clear; clc;

% From mecOfMaterials_shaft2

% sigma_max = ?
% sigma_min = ?
% sigma_mean = (sigma_max + sigma_min)/2; % Mean (midrange) stess 
% sigma_amp = (sigma_max - sigma_min)/2; % Alternating stress amplitude
% 
% M_y = ?
% M_z = ?
% 
% M = sqrt(M_y^2 + M_z^2);
%%
% Neubler's Constant for Steels, Table 6-6 page 382
S_ut_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240];
a_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009];
Neublers_table = dictionary(S_ut_table, a_table);

S_ut = S_ut_table(1, 1);
a = Neublers_table(S_ut)
% q = 1/(1 + (a/r)); % equation 6.13 page 381
%%

% Estimating stress geometric concentration factors for preliminary stage
% Shoulder fillet sharp (r/d = 0.02, D/d = 1.5), lec 5 s 10, maskinelementer)
% K_tBend = 2.7;
% K_tTor = 2.2;
% K_tAxial = 3.0;
% 
% % Fatigue concentration foactor, lec 5 slide 11
% % K_fBend = 1 + q * (K_tBend - 1); % q unknown
% % K_fTor = 1 + q * (K_tTor - 1); % q unknown
% % K_fAxial = 1 + q * (K_tAxial - 1); % q unknown
% 
% % Conservative estimate for preliminary stage
% K_fBend = K_tBend;
% K_fTor = K_tTor;
% K_fAxial = K_tAxial;
% 
% % Correction factors for preliminary stage, lec 5 s 12
% % S_e = C_load * C_size * C_surf * C_temp * C_reliab * S_ePrime;
% C_load = 1; % Bending
% C_size = 1;
% % C_surf = A * S_ut^b; % A and b from table 6-3 page 369
% % C_temp = ?
% % C_reliab = ?
% 
