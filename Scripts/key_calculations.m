clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Key calculations - Shaft 1 & 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
export_import = fullfile(pwd, 'export_import');

if exist(fullfile(export_import,'deflection_shaft1.mat'), 'file')
    load(fullfile(export_import,'deflection_shaft1.mat'))
else
    error('Load deflection.m first')
end
if exist(fullfile(export_import,'deflection_shaft3.mat'), 'file')
    load(fullfile(export_import,'deflection_shaft3.mat'))
else
    error('Load deflection.m first')
end
if exist(fullfile(export_import,'loadingDiagram_shaft1.mat'), 'file')
    load(fullfile(export_import,'loadingDiagram_shaft1.mat'), 'T_M')
else
    error('Load loadingDiagrams.m first')
end
if exist(fullfile(export_import,'loadingDiagram_shaft3.mat'), 'file')
    load(fullfile(export_import,'loadingDiagram_shaft3.mat'), 'T_out')
else
    error('Load loadingDiagrams.m first')
end

%% Parameters

d_shaft_1 = d_B; % [m] % Diameter at key placement
d_shaft_3 = d_G; % [m] % Diameter at key placement

T_shaft1 = T_M; % [Nm]
T_shaft3 = T_out; % [Nm]

S_yield = 190 * 10^6;  % [Pa] medium carbon steel typical e-modulus found online
S_yield_comp = 190 * 10^6; % [Pa]

%% Calculated values

w_1 = d_shaft_1/4;   % [m]
h_1 = d_shaft_1/6;   % [m]
l_1 = 1.5*d_shaft_1; % [m]
  
w_3 = d_shaft_3/4;   % [m]
h_3 = d_shaft_3/6;   % [m]
l_3 = 1.5*d_shaft_3; % [m]


%% Check for shear failure - one key

r_shaft_1 = d_shaft_1/2; % [m]
tau_key_1 = (T_shaft1/r_shaft_1)/(w_1*l_1); % [Pa]
n_shear_1 = (0.577*S_yield)/tau_key_1; % [-]


r_shaft_3 = d_shaft_3/2; % [m]
tau_key_3 = (T_shaft3/r_shaft_3)/(w_3*l_3); % [Pa]
n_shear_3 = (0.577*S_yield)/tau_key_3; % [-]

%% Check for shear failure - two keys

r_shaft_1 = d_shaft_1/2; % [m]
tau_key_1_2 = (T_shaft1/r_shaft_1)/ (2*(w_1*l_1)); % [Pa]
n_shear_1_2 = (0.577*S_yield)/tau_key_1_2; % [-]


r_shaft_3 = d_shaft_3/2; % [m]
tau_key_3_2 = (T_shaft3/r_shaft_3)/  (2*(w_3*l_3)); % [Pa]
n_shear_3_2 = (0.577*S_yield)/tau_key_3_2; % [-]


%% Check for compression failure - one key

sigma_key_1 = (T_shaft1/r_shaft_1)   /      ((h_1/2)*l_1); % [Pa]
n_compression_1 = (S_yield_comp)   /       sigma_key_1; % [-]

sigma_key_3 = (T_shaft3/r_shaft_3)      /   ((h_3/2)*l_3); % [Pa]
n_compression_3 = (S_yield_comp)      /    sigma_key_3; % [-]

%% Check for compression failure - two keys

sigma_key_1_2 = (T_shaft1/r_shaft_1)   /      (2*((h_1/2)*l_1)); % [Pa]
n_compression_1_2 = (S_yield_comp)   /       sigma_key_1_2; % [-]

sigma_key_3_2 = (T_shaft3/r_shaft_3)      /   (2*((h_3/2)*l_3)); % [Pa]
n_compression_3_2 = (S_yield_comp)      /    sigma_key_3_2; % [-]

%% Print results
fprintf('\n--- Shear Failure Check (One Key) ---\n');
fprintf('Shear stress key 1 (tau_key_1): %.6f [Pa]\n', tau_key_1);
fprintf('Shear safety factor key 1 (n_shear_1): %.2f\n', n_shear_1);
fprintf('Shear stress key 3 (tau_key_3): %.2f [MPa]\n', tau_key_3/ 1e6);
fprintf('Shear safety factor key 3 (n_shear_3): %.2f\n', n_shear_3);

fprintf('\n--- Shear Failure Check (Two Keys) ---\n');
fprintf('Shear stress key 1 (tau_key_1_2): %.6f [Pa]\n', tau_key_1_2);
fprintf('Shear safety factor key 1 (n_shear_1_2): %.2f\n', n_shear_1_2);
fprintf('Shear stress key 3 (tau_key_3_2): %.2f [MPa]\n', tau_key_3_2/1e6);
fprintf('Shear safety factor key 3 (n_shear_3_2): %.2f\n', n_shear_3_2);

fprintf('\n--- Compression Failure Check (One Key) ---\n');
fprintf('Compression stress key 1 (sigma_key_1): %.6f [Pa]\n', sigma_key_1);
fprintf('Compression safety factor key 1 (n_compression_1): %.2f\n', n_compression_1);
fprintf('Compression stress key 3 (sigma_key_3): %.2f [MPa]\n', sigma_key_3/1e6);
fprintf('Compression safety factor key 3 (n_compression_3): %.2f\n', n_compression_3);

fprintf('\n--- Compression Failure Check (Two Keys) ---\n');
fprintf('Compression stress key 1 (sigma_key_1_2): %.6f [Pa]\n', sigma_key_1_2);
fprintf('Compression safety factor key 1 (n_compression_1_2): %.2f\n', n_compression_1_2);
fprintf('Compression stress key 3 (sigma_key_3_2): %.2f [MPa]\n', sigma_key_3_2/1e6);
fprintf('Compression safety factor key 3 (n_compression_3_2): %.2f\n', n_compression_3_2);
