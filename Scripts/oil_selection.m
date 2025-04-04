clc; close all; clear;

export_import = fullfile(pwd, 'export_import');

if exist(fullfile(export_import,'loadingDiagram_common.mat'), 'file')
    load(fullfile(export_import,'loadingDiagram_common.mat'), ...
    "F_t1","F_t2","F_t3","F_t4","i_s1","i_s2")
end

if exist(fullfile(export_import,'gear_sizes.mat'), 'file')
    load(fullfile(export_import,'gear_sizes.mat'), ...
    "b_s1", "b_s2", "d_g3", "d_g4", "d_g1", "d_g2", ...
    "v_p_g1","v_p_g2", "v_p_g3", "v_p_g4",...
    "n_1", "n_2", "n_3","n_4")
end

if exist(fullfile(export_import,'bearings.mat'), 'file')
    load(fullfile(export_import,'bearings.mat'), ...
    "dm_B","dm_C","dm_D","dm_E","dm_F","dm_G")
end

%% Gear oil

% Gear one
loadSpeed_factor_g1 = ks(F_t1, b_s1, d_g1,i_s1) / v_p_g1; % [N·s]/[mm²·m]
% Gear two
loadSpeed_factor_g2 = ks(F_t2, b_s1, d_g2,i_s1) / v_p_g2; % [N·s]/[mm²·m]
% Gear three
loadSpeed_factor_g3 = ks(F_t3, b_s2, d_g3,i_s2) / v_p_g3; % [N·s]/[mm²·m]
% Gear four
loadSpeed_factor_g4 = ks(F_t4, b_s2, d_g4,i_s2) / v_p_g4; % [N·s]/[mm²·m]

format compact
fprintf('Load and speed factor for gear 1: %.1f [N·s]/[mm²·m]\n', loadSpeed_factor_g1);
fprintf('Load and speed factor for gear 2: %.1f [N·s]/[mm²·m]\n', loadSpeed_factor_g2);
fprintf('Load and speed factor for gear 3: %.1f [N·s]/[mm²·m]\n', loadSpeed_factor_g3);
fprintf('Load and speed factor for gear 4: %.1f [N·s]/[mm²·m]\n', loadSpeed_factor_g4);
fprintf('\n');

%% Bearing Oil

% Bearing B
fprintf('Mean diameter for bearing B: %.1f [mm]\n', dm_B);
fprintf('RPM of shaft 1: %.1f [rpm]\n', n_1);
fprintf('\n');

% Bearing C
fprintf('Mean diameter for bearing C: %.1f [mm]\n', dm_C);
fprintf('RPM of shaft 1: %.1f [rpm]\n', n_1);
fprintf('\n');

% Bearing D
fprintf('Mean diameter for bearing D: %.1f [mm]\n', dm_D);
fprintf('RPM of shaft 2: %.1f [rpm]\n', n_2);
fprintf('\n');

% Bearing E
fprintf('Mean diameter for bearing E: %.1f [mm]\n', dm_E);
fprintf('RPM of shaft 2: %.1f [rpm]\n', n_2);
fprintf('\n');

% Bearing F
fprintf('Mean diameter for bearing F: %.1f [mm]\n', dm_F);
fprintf('RPM of shaft 3: %.1f [rpm]\n', n_4);
fprintf('\n');

% Bearing G
fprintf('Mean diameter for bearing G: %.1f [mm]\n', dm_G);
fprintf('RPM of shaft 3: %.1f [rpm]\n', n_4);
fprintf('\n');


%% Functions 
function ks = ks(F_t, b, d1, U)
    % Beregner k_s til Belastning-hastighetsfaktorer
    % Input:
    %   Ft - Tangential Gear Force [N]
    %   b  - Gear Width [mm]
    %   d1 - Pitch Circle Diameter [mm]
    %   U  - Gear Ratio

    % Output:
    %   ks - [N/mm^2]

    ks = (F_t / (b * d1)) * ((U + 1) / U) * 3;
end
