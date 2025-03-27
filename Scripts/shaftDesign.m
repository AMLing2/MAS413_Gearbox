close all; clc;


% TO DO:
% - Complete check if stress due to shear can be necglected
% - Interpolation for D/d other than 1.2


% Common input parameters (for all shafts)
n_f = 2; % Safety factor
material = 355; % (355, 4140)
load_type = "Complex axial"; % ("Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial");
surface_finish = "Machined"; % ("Ground" "Machined" "Hot-rolled" "As-forged") Other types: Machine Design pg 368, fig 6-26
reliability = 95; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)
operating_temperature = 70; % [celsius] defined by Kjell (only significant if > 450)

% Cross sections lists
% cs_ = [diameter (mm), notch fillet radius (mm), D_d, P_x (N), M_z(Nmm), M_y (Nmm), T (Nmm),...
%        1|0 (1 for keyseat 0 for shoulder-fillet), K_t_keyseat]


%%%%%%%%%%%% Shaft 1 %%%%%%%%%%%%
r_keyseat1 = 0.25; % [mm]
K_t_keyseat1 = 3.75; % Keyseat stress concentration factor % Machine Design fig 10-16 pg 615
r_fillet1 = 2; % [mm]
D_d_1 = 1.2; % [-] ASK ABOUT THIS - RKH

% Import
load('loadingDiagram_shaft1.mat', 'cs_A', 'cs_0', 'cs_1', 'cs_2', ...
    'd_B', 'd_S1', 'd_12', 'd_C')

% [m] to [mm]
d_B  = d_B*1000;  % [mm]
d_S1 = d_S1*1000; % [mm]
d_12 = d_12*1000; % [mm]
d_C  = d_C*1000;  % [mm]

cs_A = [cs_A d_B r_keyseat1 D_d_1 1 K_t_keyseat1];
cs_0 = [cs_0 d_B r_fillet1  D_d_1 0 0];
cs_1 = [cs_1 d_S1 r_fillet1 D_d_1 0 0];
cs_2 = [cs_2 d_C r_fillet1  D_d_1 0 0];

%%%%%%%%%%%% Shaft 2 %%%%%%%%%%%%
r_fillet2 = 2; % [mm]
D_d_2 = 1.2; % [-] ASK ABOUT THIS - RKH

% Import
load('loadingDiagram_shaft2.mat', 'cs_3', 'cs_4', 'cs_5', 'cs_6', ...
    'd_E', 'd_S22', 'd_45', 'd_S21', 'd_D')

% [m] to [mm]
d_E = d_E*1000;     % [mm]
d_S22 = d_S22*1000; % [mm]
d_45 = d_45*1000;   % [mm]
d_S21 = d_S21*1000; % [mm]
d_D = d_D*1000;     % [mm]

cs_6 = [cs_6 d_E    r_fillet2 D_d_2 0 0];
cs_5 = [cs_5 d_S22  r_fillet2 D_d_2 0 0];
cs_4 = [cs_4 d_S21  r_fillet2 D_d_2 0 0];
cs_3 = [cs_3 d_D    r_fillet2 D_d_2 0 0];

%%%%%%%%%%%% Shaft 3 %%%%%%%%%%%%
r_keyseat3 = 0.25; % [mm]
K_t_keyseat3 = 3.75; % Keyseat stress concentration factor % Machine Design fig 10-16 pg 615
r_fillet3 = 2; % [mm]
D_d_3 = 1.2; % !! ASK ABOUT THIS !!

% Import
load('loadingDiagram_shaft3.mat', 'cs_7', 'cs_8', 'cs_9', 'cs_H', ...
    'd_F', 'd_78', 'd_S3', 'd_G')

% [m] to [mm]
d_F = d_F*1000;   % [mm]
d_78 = d_78*1000; % [mm]
d_S3 = d_S3*1000; % [mm]
d_G = d_G*1000;   % [mm]

cs_7 = [cs_7 d_F  r_fillet3  D_d_3 0 0];
cs_8 = [cs_8 d_78 r_fillet3  D_d_3 0 0];
cs_9 = [cs_9 d_S3 r_fillet3  D_d_3 0 0];
cs_H = [cs_H d_G  r_keyseat3 D_d_3 1 K_t_keyseat3];


%%%%%%%%%%%%%%%% Select cross section to evaluate %%%%%%%%%%%%%%%%
cs = cs_A;
first_iteration = "n";  % ("y" / "n") First iteration for diameter equation (limited geometry data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lists for for loop
shaft1_AllCS = [cs_A; cs_0; cs_1; cs_2];
shaft1_names = {'Cross section A', 'Cross section 0', 'Cross section 1', 'Cross section 2'};
shaft2_AllCS = [cs_6; cs_5; cs_4; cs_3];
shaft2_names = {'Cross section 6', 'Cross section 5', 'Cross section 4', 'Cross section 3'};
shaft3_AllCS = [cs_7; cs_8; cs_9; cs_H];
shaft3_names = {'Cross section 7', 'Cross section 8', 'Cross section 9', 'Cross section H'};


%% Shaft 1 loop

shaft1_results = zeros(size(shaft1_AllCS, 1), 5);  % 4 x 5
for i = 1:size(shaft1_AllCS, 1)
    fprintf('\n------ %s ------\n', shaft1_names{i});
    cs = shaft1_AllCS(i, :);
    run("fatigue.m")
    shaft1_results(i, :) = [S_e, sigma_vm_mean, sigma_vm_amp, n_y, n_f];
end

modifiedGoodman(S_y, S_yc, S_ut, shaft1_results, 'Shaft 1')

%% Shaft 2 loop

shaft2_results = zeros(size(shaft2_AllCS, 1), 5);  % 4 x 5
for i = 1:size(shaft2_AllCS, 1)
    fprintf('\n------ %s ------\n', shaft2_names{i});
    cs = shaft2_AllCS(i, :);
    run("fatigue.m")
    shaft2_results(i, :) = [S_e, sigma_vm_mean, sigma_vm_amp, n_y, n_f];
end

modifiedGoodman(S_y, S_yc, S_ut, shaft2_results, 'Shaft 2')

%% Shaft 3 loop

shaft3_results = zeros(size(shaft3_AllCS, 1), 5);  % 4 x 5
for i = 1:size(shaft3_AllCS, 1)
    fprintf('\n------ %s ------\n', shaft3_names{i});
    cs = shaft3_AllCS(i, :);
    run("fatigue.m")
    shaft3_results(i, :) = [S_e, sigma_vm_mean, sigma_vm_amp, n_y, n_f];
end

modifiedGoodman(S_y, S_yc, S_ut, shaft3_results, 'Shaft 3')