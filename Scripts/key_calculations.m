clc; close all; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Key calculations - Shaft 1 & 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('deflection_shaft1.mat')
load('deflection_shaft3.mat')

load('loadingDiagrams_shaft1.mat', "T_M")
load('loadingDiagrams_shaft3.mat', "T_out")

%% Parameters

d_shaft_1 = d_B; %[m] % Diameter at key placement
d_shaft_3 = d_G; %[m] % Diameter at key placement

T_shaft1 = T_M; %[Nm]
T_shaft3 = T_out; %[Nm]

S_yield = 190 * 10^9;  % [Pa] medium carbon steel typical e-modulus found online

%% Calculated values

r_shaft_1 = d_shaft_1/2; %[m]
r_shaft_3 = d_shaft_3/2; %[m]

w_1 = d_shaft_1/4;  %[Nm]
h_1 = d_shaft_1/6;  %[Nm]
l_1 = 1.5*d_shaft_1;%[Nm]

w_3 = d_shaft_3/4;  %[Nm]
h_3 = d_shaft_3/6;  %[Nm]
l_3 = 1.5*d_shaft_3;%[Nm]


%% Check for shear failure

tau_key_1 = (T_shaft1/r_shaft_1)/(w_1*l_1); %[Pa]
n_shear_1 = (0.577*S_yield)/tau_key_1; %[-]

tau_key_3 = (T_shaft3/r_shaft_3)/(w_3*l_3); %[Pa]
n_shear_3 = (0.577*S_yield)/tau_key_3; %[-]


%% Check for compression failure

sigma_key_1 = (T_shaft1/r_shaft_1)/((h_1/2)*l_1); %[Pa]
n_compression_1 = (0.577*S_yield)/sigma_key_1; %[-]

sigma_key_3 = (T_shaft3/r_shaft_3)/((h_3/2)*l_3); %[Pa]
n_compression_3 = (0.577*S_yield)/sigma_key_3; %[-]

%% Print results
fprintf('\n--- Shear Failure Check ---\n');
fprintf('Shear stress key 1 (tau_key_1): %.2f Pa\n', tau_key_1);
fprintf('Shear safety factor key 1 (n_shear_1): %.2f\n', n_shear_1);
fprintf('Shear stress key 3 (tau_key_3): %.2f Pa\n', tau_key_3);
fprintf('Shear safety factor key 3 (n_shear_3): %.2f\n', n_shear_3);

fprintf('\n--- Compression Failure Check ---\n');
fprintf('Compression stress key 1 (sigma_key_1): %.2f Pa\n', sigma_key_1);
fprintf('Compression safety factor key 1 (n_compression_1): %.2f\n', n_compression_1);
fprintf('Compression stress key 3 (sigma_key_3): %.2f Pa\n', sigma_key_3);
fprintf('Compression safety factor key 3 (n_compression_3): %.2f\n', n_compression_3);