clc; clear; close all;
export_import = fullfile(pwd, 'export_import');

% Common input parameters (for all shafts)
n_f = 2; % Safety factor
S355J2 = [315, 470, 210*1e9, 0.3]; % Material data [S_y, S_ut, Youngs module (Pa), Poisson´s ratio]
material = 355; % (355, 4140)
load_type = "Complex axial"; % ("Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial");
surface_finish = "Machined"; % ("Ground" "Machined" "Hot-rolled" "As-forged") Other types: Machine Design pg 368, fig 6-26
reliability = 99; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)
operating_temperature = 70; % [celsius] defined by Kjell (only significant if > 450)
first_iteration = false;  % (true / false) First iteration for diameter equation (limited geometry data)

% % Shaft Diameters (Old values)
% d_C   = 0.010; % [m]
% d_12  = 0.015; % [m]
% d_B   = 0.011; % [m]
% d_S1  = 0.013; % [m]
% d_D = 0.01;    % [m]
% d_E = 0.015;   % [m]
% d_S21 = 0.02;  % [m]
% d_S22 = 0.02;  % [m]
% d_45 = 0.02;   % [m]
% d_F   = 0.020; % [m]
% d_78  = 0.030; % [m]
% d_G = 0.027;   % [m]
% d_S3 = 0.025;  % [m]


%%%%%%%%%%%% Shaft 1 %%%%%%%%%%%%
r_keyseat1 = 0.25; % [mm] Keyseat fillet radius
K_t_keyseat1 = 3.75; % Keyseat stress concentration factor % Machine Design fig 10-16 pg 615
r_fillet1 = 0.4; % [mm] Shoulder fillet radius

% Diameters
d_B  = 25.025; % [mm] 
d_S1 = 29;     % [mm]
d_12 = 44;     % [mm]
d_C  = 35.035; % [mm]

%%%%%%%%%%%% Shaft 2 %%%%%%%%%%%%
r_fillet2 = 1; % [mm] Shoulder fillet radius

% Diameters
d_S22 = 54;     % [mm]
d_E   = 45.045; % [mm]
d_45  = 60;     % [mm]
d_S21 = 52;     % [mm]
d_D   = 40.04;  % [mm]

%%%%%%%%%%%% Shaft 3 %%%%%%%%%%%%
r_keyseat3 = 0.6; % [mm] Keyseat fillet radius
K_t_keyseat3 = 3.8; % Keyseat stress concentration factor % Machine Design fig 10-16 pg 615
r_fillet3 = 0.5; % [mm] Shoulder fillet radius

% Diameters
d_G  = 70.07; % [mm]
d_S3 = 76;     % [mm]
d_78 = 82;     % [mm]
d_F  = 65.065;  % [mm]

if exist(fullfile(export_import, "shrinkFit_diameters.mat"), 'file')
    load(fullfile(export_import, 'shrinkFit_diameters.mat'), ...
        'd_S1', 'd_B', 'd_C', 'd_12', 'd_S22', 'd_E', 'd_D', ...
        'd_S21', 'd_45', 'd_78', 'd_F', 'd_G', 'd_S3')
else
    fprintf('Using diameters from shaftDesign\n')
end

% Import
% Cross sections lists
% cs_ = [P (N), T (Nmm), M (Nmm), V_y (N), V_z (N), shaft diameter,...
%        boolean (true -> keyseat, false -> shoulder-fillet),...
%        notch radius, (keyseat stress concentration | D/d)];

if exist(fullfile("export_import","loadingDiagram_shaft1.mat"),'file')
    load(fullfile("export_import", 'loadingDiagram_shaft1.mat'), 'cs_A', 'cs_0', 'cs_1', 'cs_2')
    load(fullfile("export_import", 'loadingDiagram_shaft2.mat'), 'cs_3', 'cs_4', 'cs_5', 'cs_6','E_Fa','D_Fa','L_45')
    load(fullfile("export_import", 'loadingDiagram_shaft3.mat'), 'cs_7', 'cs_8', 'cs_9', 'cs_H')
else
    error("Run loadingDiagrams.m first")
end

cs_A = [cs_A d_B   true  r_keyseat1 K_t_keyseat1];
cs_0 = [cs_0 d_B   false r_fillet1  d_S1/d_B];
cs_1 = [cs_1 d_S1  false r_fillet1  d_12/d_S1];
cs_2 = [cs_2 d_C   false r_fillet1  d_12/d_C];
cs_6 = [cs_6 d_E   false r_fillet2  d_S22/d_E];
cs_5 = [cs_5 d_S22 false r_fillet2  d_45/d_S22];
cs_4 = [cs_4 d_S21 false r_fillet2  d_45/d_S21];
cs_3 = [cs_3 d_D   false r_fillet2  d_S21/d_D];
cs_7 = [cs_7 d_F   false r_fillet3  d_78/d_F];
cs_8 = [cs_8 d_78  false r_fillet3  d_78/d_S3];
cs_9 = [cs_9 d_S3  false r_fillet3  d_S3/d_G];
cs_H = [cs_H d_G   true  r_keyseat3 K_t_keyseat3];

% cs_A  = [cs_A d_B   true  r_keyseat1 K_t_keyseat1];
% cs_0l = [cs_0 d_B   false r_fillet1  d_S1/d_B];
% cs_0r = [cs_0 d_B   false r_fillet1  d_S1/d_B];
% cs_1l = [cs_1 d_S1  false r_fillet1  d_12/d_S1];
% cs_1r = [cs_1 d_S1  false r_fillet1  d_12/d_S1];
% cs_2l = [cs_2 d_C   false r_fillet1  d_12/d_C];
% cs_2r = [cs_2 d_C   false r_fillet1  d_12/d_C];
% cs_6l = [cs_6 d_E   false r_fillet2  d_S22/d_E];
% cs_6r = [cs_6 d_E   false r_fillet2  d_S22/d_E];
% cs_5l = [cs_5 d_S22 false r_fillet2  d_45/d_S22];
% cs_5r = [cs_5 d_S22 false r_fillet2  d_45/d_S22];
% cs_4l = [cs_4 d_S21 false r_fillet2  d_45/d_S21];
% cs_4r = [cs_4 d_S21 false r_fillet2  d_45/d_S21];
% cs_3l = [cs_3 d_D   false r_fillet2  d_S21/d_D];
% cs_3r = [cs_3 d_D   false r_fillet2  d_S21/d_D];
% cs_7l = [cs_7 d_F   false r_fillet3  d_78/d_F];
% cs_7r = [cs_7 d_F   false r_fillet3  d_78/d_F];
% cs_8l = [cs_8 d_78  false r_fillet3  d_78/d_S3];
% cs_8r = [cs_8 d_78  false r_fillet3  d_78/d_S3];
% cs_9l = [cs_9 d_S3  false r_fillet3  d_S3/d_G];
% cs_9r = [cs_9 d_S3  false r_fillet3  d_S3/d_G];
% cs_H  = [cs_H d_G   true  r_keyseat3 K_t_keyseat3];

% Lists for for loop
shaft1_AllCS = [cs_A; cs_0; cs_1; cs_2];
shaft2_AllCS = [cs_3; cs_4; cs_5; cs_6];
shaft3_AllCS = [cs_7; cs_8; cs_9; cs_H];
shaft1_names = {'Cross section A', 'Cross section 0', 'Cross section 1', 'Cross section 2'};
shaft2_names = {'Cross section 3', 'Cross section 4', 'Cross section 5', 'Cross section 6'};
shaft3_names = {'Cross section 7', 'Cross section 8', 'Cross section 9', 'Cross section H'};

% shaft1_AllCS = [cs_A; cs_0l; cs_0r; cs_1l; cs_1r; cs_2l; cs_2r];
% shaft2_AllCS = [cs_3l; cs_3r; cs_4l; cs_4r; cs_5l; cs_5r; cs_6l; cs_6r];
% shaft3_AllCS = [cs_7l; cs_7r; cs_8l; cs_8r; cs_9l;  cs_9r; cs_H];
% shaft1_names = {'Cross section A', 'Cross section 0', 'Cross section 1', 'Cross section 2'};
% shaft2_names = {'Cross section 3l', 'Cross section 3r', 'Cross section 4l', 'Cross section 4r'
%                 'Cross section 5l', 'Cross section 5r', 'Cross section 6l', 'Cross section 6r'};
% shaft3_names = {'Cross section 7', 'Cross section 8', 'Cross section 9', 'Cross section H'};


% %% Shaft 1 loop
% 
% shaft1_results = zeros(size(shaft1_AllCS, 1), 5);
% for i = 1:size(shaft1_AllCS, 1)
%     fprintf('\n------ %s ------\n', shaft1_names{i});
%     cs = shaft1_AllCS(i, :);
%     run("fatigue.m")
%     shaft1_results(i, :) = [S_e, sigma_e_m, sigma_e_a, n_y, n_f];
% end
% 
% modifiedGoodman(S_y, S_ut, shaft1_results, 'Shaft 1')

%% Shaft 2 loop

shaft2_results = zeros(size(shaft2_AllCS, 1), 5);
for i = 1:size(shaft2_AllCS, 1)
    fprintf('\n\n------ %s ------\n', shaft2_names{i});
    cs = shaft2_AllCS(i, :);
    run("fatigue.m")
    shaft2_results(i, :) = [S_e, sigma_e_m, sigma_e_a, n_y, n_f];
end

modifiedGoodman(S_y, S_ut, shaft2_results, 'Shaft 2')

% %% Shaft 3 loop
% 
% shaft3_results = zeros(size(shaft3_AllCS, 1), 5);
% for i = 1:size(shaft3_AllCS, 1)
%     fprintf('\n\n------ %s ------\n', shaft3_names{i});
%     cs = shaft3_AllCS(i, :);
%     run("fatigue.m")
%     shaft3_results(i, :) = [S_e, sigma_e_m, sigma_e_a, n_y, n_f];
% end
% 
% modifiedGoodman(S_y, S_ut, shaft3_results, 'Shaft 3')
% load(fullfile(export_import, "loadingDiagram_common.mat"),'F_a_remaining')
% 
% n_f = 3; % reset n_f
% shoulder_length_min = shoulder_length(F_a_remaining,max(d_S22,d_S21),...
%     d_45-r_fillet2*1.5*2,n_f,S_y); % calculate gear chamfer not acting on the bottom of the shoulder
% if shoulder_length_min < L_45*1e3
%     fprintf("\nL_45 shoulder length good\n")
% else
%     fprintf("\nL_45 shoulder length needs to be increased to above %f [mm]\n",shoulder_length_min)
% end

% Export diameters
save(fullfile("export_import","shaftDesign.mat"), 'd_S1', 'd_B',...
'd_C', 'd_12', 'd_S22', 'd_E', 'd_D', 'd_S21', 'd_45', 'd_78', 'd_F', ...
'd_G', 'd_S3', 'r_fillet1','r_fillet2','r_fillet3', 'E','V_shaft')