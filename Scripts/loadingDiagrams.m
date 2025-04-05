clc; clear; close all;
export_import = fullfile(pwd, 'export_import');

disablePlotting = false;
disablePlot1 = true;
disablePlot2 = true;
disablePlot3 = false;

%% Constants

% Chosen Parameters
L_12 = 5e-3; % [m]
L_45 = 5e-3; % [m]
L_78 = 5e-3; % [m]
L_AB = 0.05; % [m]
L_GH = 100e-3; % [m]
    % Stage efficiencies
% eta = 0.96; % [-] "finely worked teeth & good lubrication"
eta = 1.00; % [-] Ideal Stages

% Common Plotting Constants
colFill = [0.7765 0.9176 0.9843];
resolution = 100;
wPlot = 22;
hPlot = 16;
fSize = 16;

% Given information
n_1 = 1450; % [RPM]
P_1 = 12.5e3; % [W]
i_tot_og = 17.3;
alpha = 20; % [degrees] Helix Angle
beta = 15;  % [degrees] Pressure Angle

% Gravity Constant
g = 9.81; % [m/s^2]

%% Import

% Import from Bearings
if exist(fullfile(export_import,'bearings.mat'), 'file')
    load(fullfile(export_import,'bearings.mat'), ...
        'b_B', 'b_C', 'b_E', 'b_D', 'b_F', 'b_G') % [mm]
    b_F = b_F / 1000; % [m]
    b_G = b_G / 1000; % [m]
    b_E = b_E / 1000; % [m]
    b_D = b_D / 1000; % [m]
    b_B = b_B / 1000; % [m]
    b_C = b_C / 1000; % [m]
else
    % Bearing Widths
    b_B = 30e-3; % [m]
    b_C = b_B;   % [m]
    b_E = 30e-3; % [m]
    b_D = b_E;   % [m]
    b_F = 30e-3; % [m]
    b_G = b_F;   % [m]
end

% Import from Gear Sizing
if exist(fullfile(export_import, 'gear_sizes.mat'), 'file')
    load(fullfile(export_import,'gear_sizes.mat'), ...
        'd_g1', 'd_g2', 'd_g3', 'd_g4', ...
        'b_s1', 'b_s2', 'i_tot', 'i_s1','i_s2', ...
        'mass_g1', 'mass_g2', 'mass_g3', 'mass_g4') % [mm], [MPa], [kg]
    r_G1 = d_g1/2 * 1e-3; % [m]
    r_G2 = d_g2/2 * 1e-3; % [m]
    r_G3 = d_g3/2 * 1e-3; % [m]
    r_G4 = d_g4/2 * 1e-3; % [m]
    b_s1 = b_s1 * 1e-3; % [m]
    b_s2 = b_s2 * 1e-3; % [m]
else
    error('Start by running gear_sizing.m')
end

%% Common Calculations

% General
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
n_out = (n_1/i_tot); % [RPM]
omega_out = n_out * 2*pi / 60; % [rad/sec]
eta_tot = eta^2; % [-] Squared because there are two stages
P_out = P_1*eta_tot; % [W]
T_M   = P_1/omega_1; % [Nm]
T_out = P_out/omega_out; %[Nm]

% Gear 1
F_t1 = T_M / r_G1; % [N]
F_a1 = F_t1 * tand(beta); % [N]
F_r1 = F_t1 * tand(alpha)/cosd(beta); % [N]
F_G1 = mass_g1*g; % [N]
% Gear 2
F_t2 = T_M/r_G1; % [N]
F_a2 = F_t2 * tand(beta); % [N]
F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]
F_G2 = mass_g2*g;  % [N]
% Gear 3
F_t3 = T_out / r_G4; % [N]
F_a3 = F_t3 * tand(beta); % [N]
F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]
F_G3 = mass_g3*g; % [N]
% Gear 4
F_t4 = T_out / r_G4; % [N]
F_a4 = F_t4 * tand(beta); % [N]
F_r4 = F_t4 * tand(alpha)/cosd(beta); % [N]
F_G4 = mass_g4*g; % [N]

% Remaining Axial force of shaft 2
F_a_remaining = abs(F_a3 - F_a2)

%% Export Common for all Shafts
save(fullfile(export_import, "loadingDiagram_common.mat"))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAS413 Project: Loading Diagrams - Shaft 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lengths
L_BC  = b_B/2 + L_78 + b_s2 + L_45 + b_s1 + L_12 + b_C/2; % [m]
L_BG1 = b_B/2 + L_78 + b_s2 + L_45 + b_s1/2; % [m]
L_G1C = L_BC - L_BG1; % [m]
L_AG1 = L_AB + L_BG1; % [m]
L_AC = L_AG1 + L_G1C; % [m]
L_A0 = L_AB + b_B/2;  % [m]
L_A1 = L_AG1 + b_s1/2; % [m]
L_A2 = L_A1 + L_12; % [m]
% For Reaction forces @ bearings
F_By = F_t1*L_G1C/L_BC; % [N]
F_Bz = (F_r1*L_G1C - F_a1*r_G1)/L_BC; % [N]
F_Bg1 = (F_G1*L_G1C)/L_BC;


%% XY - Plane

% Figure setup
figHandle = 1;
xPos = 10;
yPos = 3;

% Initialize XY Plots
xy_x = [];
xy_P = [];
xy_V = [];
xy_M = [];
xy_T = [];

% 0 < x < L_AB
x = linspace(0, L_AB, resolution);
x(size(x))
xy_x = [xy_x, x]; % [m]
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, zeros(size(x))]; % [N]
xy_M = [xy_M, zeros(size(x))]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M ]; % [Nm]

% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xy_x = [xy_x, x]; % [m]
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, ones(size(x)) * F_By]; % [N]
xy_M = [xy_M, ( F_By*(x - L_AB) )]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M]; % [Nm]

% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xy_x = [xy_x, x]; % [m]
xy_P = [xy_P, ones(size(x)) * (-F_a1)]; % [N]
xy_V = [xy_V, ones(size(x)) * (F_By - F_t1)]; % [N]
xy_M = [xy_M, ( F_By*(x - L_AB) - F_t1*(x - L_AG1) )]; % [Nm]
xy_T = [xy_T, ones(size(x)) * (T_M - F_t1*r_G1)]; % [Nm]

if (~disablePlotting) && (~disablePlot1)
    figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 1: XY - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_y(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_z(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')

    subplot(2,2,1)
    plotLD(xy_x,xy_P, colFill)
    subplot(2,2,2)
    plotLD(xy_x,xy_V, colFill)
    subplot(2,2,3)
    plotLD(xy_x,xy_M, colFill)
    subplot(2,2,4)
    plotLD(xy_x,xy_T, colFill)
end

%% XZ - Plane

% Figure setup
figHandle = 2;
xPos = 10;
yPos = 3;

% Initialize XZ Plots
xz_x = [];
xz_P = [];
xz_V = [];
xz_M = [];
xz_T = [];
xz_Mg = [];

% 0 < x < L_AB
x = linspace(0, L_AB, resolution); % L_AB blir duplikert i xz_[n] listene
xz_x = [xz_x, x]; % [m]
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, zeros(size(x))]; % [N]
xz_M = [xz_M, zeros(size(x))]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]
xz_Mg =[xz_Mg, zeros(size(x))]; % [Nm]

% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xz_x = [xz_x, x]; % [m]
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, ones(size(x)) * (-F_Bz)]; % [N]
xz_M = [xz_M, ( - F_Bz*(x - L_AB) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]
xz_Mg =[xz_Mg, +( F_Bg1 * (x-L_AB)) ]; % [Nm]

% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xz_x = [xz_x, x]; % [m]
xz_P = [xz_P, ones(size(x)) * (-F_a1)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r1 - F_Bz)]; % [N]
xz_M = [xz_M, ( F_r1*(x - L_AG1) - F_Bz*(x - L_AB) - F_a1*(r_G1) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (T_M - F_t1*r_G1)]; % [Nm]
xz_Mg =[xz_Mg, +( F_Bg1 * (x - L_AB) ) - ( F_G1 * (x - L_AG1) ) ]; % [Nm]

if (~disablePlotting) && (~disablePlot1)
    figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 1: XZ - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_z(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_y(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')
    
    subplot(2,2,1)
    plotLD(xz_x,xz_P,colFill)
    subplot(2,2,2)
    plotLD(xz_x,xz_V,colFill)
    subplot(2,2,3)
    plotLD(xz_x,xz_M,colFill)
    subplot(2,2,4)
    plotLD(xz_x,xz_T,colFill)
end

%% Loading Diagram: Combined Bending Moment

M = sqrt(xz_M.^2 + xy_M.^2); % [Nm]
[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);

% Figure setup
figHandle = 3;
wPlotM = wPlot;
hPlotM = hPlot/2;

if (~disablePlotting) && (~disablePlot1)
    figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos+hPlotM/2 wPlotM hPlotM]);
    plotLD(xy_x, M, colFill)
    title('\textbf{Shaft 1: Combined Bending Moment}', ...
        'interpreter', 'latex', 'FontSize',fSize)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    xlim([0 L_AC])
    
    dashLineV(L, figHandle, 2, 2)
end


%% Export
idx = 2;
%%%% For Fatigue %%%%
% Cross Section A
[~, cs_A_idx] = closest(xy_x, 0);
cs_A_P = xy_P(cs_A_idx); % [N]
cs_A_T = xy_T(cs_A_idx)*1e3; % [Nmm]
cs_A_M = M(cs_A_idx)*1e3; % [Nmm]
cs_A_Vy = xy_V(cs_A_idx); % [N]
cs_A_Vz = xz_V(cs_A_idx); % [N]
cs_A = [cs_A_P cs_A_T cs_A_M cs_A_Vy cs_A_Vz];
% Cross Section 0
[~, cs_0L_idx] = closest(xy_x, L_A0);
cs_0L_idx = cs_0L_idx - idx;
cs_0L_P = xy_P(cs_0L_idx); % [N]
cs_0L_T = xy_T(cs_0L_idx)*1e3; % [Nmm]
cs_0L_M = M(cs_0L_idx)*1e3; % [Nmm]
cs_0L_Vy = xy_V(cs_0L_idx); % [N]
cs_0L_Vz = xz_V(cs_0L_idx); % [N]
cs_0L = [cs_0L_P cs_0L_T cs_0L_M cs_0L_Vy cs_0L_Vz];
[~, cs_0R_idx] = closest(xy_x, L_A0);
cs_0R_idx = cs_0R_idx + idx;
cs_0R_P = xy_P(cs_0R_idx); % [N]
cs_0R_T = xy_T(cs_0R_idx)*1e3; % [Nmm]
cs_0R_M = M(cs_0R_idx)*1e3; % [Nmm]
cs_0R_Vy = xy_V(cs_0R_idx); % [N]
cs_0R_Vz = xz_V(cs_0R_idx); % [N]
cs_0R = [cs_0R_P cs_0R_T cs_0R_M cs_0R_Vy cs_0R_Vz];
% Cross Section 1
[~, cs_1L_idx] = closest(xy_x, L_A1);
cs_1L_idx = cs_1L_idx - idx;
cs_1L_P = xy_P(cs_1L_idx); % [N]
cs_1L_T = xy_T(cs_1L_idx)*1e3; % [Nmm]
cs_1L_M = M(cs_1L_idx)*1e3; % [Nmm]
cs_1L_Vy = xy_V(cs_1L_idx); % [N]
cs_1L_Vz = xz_V(cs_1L_idx); % [N]
cs_1L = [cs_1L_P cs_1L_T cs_1L_M cs_1L_Vy cs_1L_Vz];
[~, cs_1R_idx] = closest(xy_x, L_A1);
cs_1R_idx = cs_1R_idx + idx;
cs_1R_P = xy_P(cs_1R_idx); % [N]
cs_1R_T = xy_T(cs_1R_idx)*1e3; % [Nmm]
cs_1R_M = M(cs_1R_idx)*1e3; % [Nmm]
cs_1R_Vy = xy_V(cs_1R_idx); % [N]
cs_1R_Vz = xz_V(cs_1R_idx); % [N]
cs_1R = [cs_1R_P cs_1R_T cs_1R_M cs_1R_Vy cs_1R_Vz];
% Cross Section 2
[~, cs_2L_idx] = closest(xy_x, L_A2);
cs_2L_idx = cs_2L_idx - idx;
cs_2L_P = xy_P(cs_2L_idx); % [N]
cs_2L_T = xy_T(cs_2L_idx)*1e3; % [Nmm]
cs_2L_M = M(cs_2L_idx)*1e3; % [Nmm]
cs_2L_Vy = xy_V(cs_2L_idx); % [N]
cs_2L_Vz = xz_V(cs_2L_idx); % [N]
cs_2L = [cs_2L_P cs_2L_T cs_2L_M cs_2L_Vy cs_2L_Vz];
[~, cs_2R_idx] = closest(xy_x, L_A2);
cs_2R_idx = cs_2R_idx + idx;
cs_2R_P = xy_P(cs_2R_idx); % [N]
cs_2R_T = xy_T(cs_2R_idx)*1e3; % [Nmm]
cs_2R_M = M(cs_2R_idx)*1e3; % [Nmm]
cs_2R_Vy = xy_V(cs_2R_idx); % [N]
cs_2R_Vz = xz_V(cs_2R_idx); % [N]
cs_2R = [cs_2R_P cs_2R_T cs_2R_M cs_2R_Vy cs_2R_Vz];
%%%% For Fatigue %%%%
%%%% For Bearings %%%%
[~, B_idx] = closest(xy_x, L_AB);
B_Fa = xy_P(B_idx); % [N] Axial Force
B_Fr = sqrt( xy_V(B_idx+2)^2 + xz_V(B_idx+2)^2 ); % [N] Radial Force
[~, C_idx] = closest(xy_x, L_AC);
C_Fa = xy_P(C_idx); % [N] Axial Force
C_Fr = sqrt( xy_V(C_idx)^2 + xz_V(C_idx)^2 ); % [N] Radial Force
%%%% For Bearings %%%%

% Export data w/o figures (https://stackoverflow.com/questions/
                            % 45560181/avoid-saving-of-graphics-in-matlab)
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save(fullfile(export_import, "loadingDiagram_shaft1.mat"), saveVars{:})

%% Length sanity check
lW = 3;

if (~disablePlotting) && (~disablePlot1)
    figHandle = 11;
    figure(figHandle)
    hold on
    plot( [0 L_AC],  [0  0], 'k', 'LineWidth', lW)
    plot( [0 L_AG1], [1  1], 'LineWidth', lW)
    plot( [0 L_AB],  [-1 -1], 'LineWidth', lW)
    plot( [L_AB (L_AB + L_BG1)], [-1 -1], 'LineWidth', lW)
    plot( [L_AG1 (L_AG1 + L_G1C)], [-1 -1], 'LineWidth', lW)
    xlabel('Length [m]')
    legend('$L_{AC}$', '$L_{AG1}$', '$L_{AB}$', '$L_{BG1}$', '$L_{G1C}$', ...
        'Location','southoutside', 'interpreter', 'latex')
    xlim( [ (-L_AC*0.1), (L_AC + L_AC*0.1) ] )
    ylim( [-5 5] )
    title('Shaft 1: One Directional Length', 'interpreter', 'latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAS413 Project: Loading Diagrams - Shaft 2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clearvars -except export_import;
load(fullfile(export_import, "loadingDiagram_common.mat"))

% Lengths
L_ED = b_E/2 + L_78 + b_s2 + L_45 + b_s1 + L_12 + b_D/2; % [m]
L_EG3 = b_E/2 + L_78 + b_s2/2; % [m]
L_EG2 = L_EG3 + b_s2/2 + L_45 + b_s1/2; % [m]
L_G3G2 = L_EG2 - L_EG3; % [m]
L_G3D = L_ED - L_EG3; % [m]
L_G2D = L_ED - L_EG2; % [m]
L_E3 = L_ED - b_D/2; % [m]
L_E4 = L_EG2 - b_s1/2; % [m]
L_E5 = L_EG3 + b_s2/2; % [m]
L_E6 = b_E/2; % [m]

% Reaction forces @ bearings
F_Dx = F_a2 - F_a3; % to the right
F_Ex = 0;
F_Ez = (F_r3*L_G3D - F_r2*L_G2D - F_a3*r_G3 - F_a2*r_G2)/L_ED;
F_Ey = (F_t3*L_G3D + F_t2*L_G2D)/L_ED;

F_EG = ( (F_G3*L_G3D) + ( F_G2*L_G2D) ) / (L_ED);

 
%% XY - Plane
 
% Figure setup
figHandle = 4;
xPos = 10;
yPos = 3;

% Initialize XY Plots
xy_x = [];
xy_P = [];
xy_V = [];
xy_M = [];
xy_T = [];

% 0 < x < L_EG3
x = linspace(0, L_EG3, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Ex)]; % [N]
xy_V = [xy_V, ones(size(x)) * (-F_Ey)]; % [N]
xy_M = [xy_M, ( - F_Ey*(x) )]; % [Nm]
xy_T = [xy_T, zeros(size(x))]; % [Nm]

% L_EG3 < x < L_EG2
x = linspace(L_EG3, L_EG2, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Ex - F_a3)]; % [N]
xy_V = [xy_V, ones(size(x)) * (-F_Ey + F_t3)]; % [N]
xy_M = [xy_M, ( - F_Ey*(x) + F_t3*(x - L_EG3) )]; % [Nm]
xy_T = [xy_T, ones(size(x)) * F_t3*r_G3]; % [Nm]

% L_EG2 < x < L_ED
x = linspace(L_EG2, L_ED, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Ex - F_a3 + F_a2)]; % [N]
xy_V = [xy_V, ones(size(x)) * (F_t2 - F_Ey + F_t3)]; % [N]
xy_M = [xy_M, ( F_t2*(x - L_EG2) - F_Ey*(x) + F_t3*(x - L_EG3) )]; % [Nm]
xy_T = [xy_T, ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2)]; % [Nm]
 
if (~disablePlotting) && (~disablePlot2)
    XYplaneFig = figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 2: XY - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_y(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_z(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')
    
    subplot(2,2,1)
    plotLD(xy_x,xy_P, colFill)
    subplot(2,2,2)
    plotLD(xy_x,xy_V, colFill)
    subplot(2,2,3)
    plotLD(xy_x,xy_M, colFill)
    subplot(2,2,4)
    plotLD(xy_x,xy_T, colFill)
end

%% XZ - Plane

% Figure setup
figHandle = 5;
xPos = 10;
yPos = 3;

% Initialize XZ Plots
xz_x = [];
xz_P = [];
xz_V = [];
xz_M = [];
xz_T = [];
xz_Mg = [];
 
% 0 < x < L_EG3
x = linspace(0, L_EG3, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x))*(-F_Ex)]; % [N]
xz_V = [xz_V, ones(size(x))*(-F_Ez)]; % [N]
xz_M = [xz_M, ( -F_Ez*(x) )]; % [Nm]
xz_T = [xz_T, zeros(size(x))]; % [Nm]
xz_Mg = [xz_Mg, F_EG*x]; % [N]

% L_EG3 < x < L_EG2
x = linspace(L_EG3, L_EG2, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x))* ( - F_Ex - F_a3)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r3 - F_Ez)]; % [N]
xz_M = [xz_M, ( -F_Ez*(x) + F_r3*(x - L_EG3) - F_a3*(r_G3) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * F_t3*r_G3]; % [Nm]
xz_Mg = [xz_Mg, F_EG*x - F_G3*(x-L_EG3)]; % [N]

% L_EG2 < x < L_ED
x = linspace(L_EG2, L_ED, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (F_a2 - F_a3 - F_Ex)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r3 - F_Ez - F_r2)]; % [N]
xz_M = [xz_M, ( -F_Ez*(x) + F_r3*(x - L_EG3) - F_r2*(x - L_EG2) - ...
                 F_a2*(r_G2) - F_a3*(r_G3) )]; % [Nm]
xz_T = [xz_T,  ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2)]; % [Nm]
xz_Mg = [xz_Mg, F_EG*x - F_G3*(x-L_EG3) - F_G2* (x-L_EG2)]; % [N]
 
if (~disablePlotting) && (~disablePlot2)
    XZplaneFig = figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 2: XZ - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_z(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_y(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')
    
    subplot(2,2,1)
    plotLD(xz_x,xz_P,colFill)
    subplot(2,2,2)
    plotLD(xz_x,xz_V,colFill)
    subplot(2,2,3)
    plotLD(xz_x,xz_M,colFill)
    subplot(2,2,4)
    plotLD(xz_x,xz_T,colFill)
end


%% Loading Diagram: Combined Bending Moment

M = sqrt(xz_M.^2 + xy_M.^2); % [Nm]
[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);

% Figure setup
figHandle = 6;
wPlotM = wPlot;
hPlotM = hPlot/2;

if (~disablePlotting) && (~disablePlot2)
    figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos+hPlotM/2 wPlotM hPlotM]);
    plotLD(xy_x, M, colFill)
    title('\textbf{Shaft 2: Combined Bending Moment}', ...
        'interpreter', 'latex', 'FontSize',fSize)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    xlim([0 L_ED])

    dashLineV(L, figHandle, 2, 2)
end

%% Export
idx = 2;
%%%% For Fatigue %%%%
% Cross Section 3
[~, cs_3L_idx] = closest(xy_x, L_E3);
cs_3L_idx = cs_3L_idx - idx;
cs_3L_P = xy_P(cs_3L_idx);
cs_3L_T = xy_T(cs_3L_idx)*1e3; % [Nmm]
cs_3L_M = M(cs_3L_idx)*1e3; % [Nmm]
cs_3L_Vy = xy_V(cs_3L_idx);
cs_3L_Vz = xz_V(cs_3L_idx);
cs_3L = [cs_3L_P cs_3L_T cs_3L_M cs_3L_Vy cs_3L_Vz];
[~, cs_3R_idx] = closest(xy_x, L_E3);
cs_3R_idx = cs_3R_idx + idx;
cs_3R_P = xy_P(cs_3R_idx);
cs_3R_T = xy_T(cs_3R_idx)*1e3; % [Nmm]
cs_3R_M = M(cs_3R_idx)*1e3; % [Nmm]
cs_3R_Vy = xy_V(cs_3R_idx);
cs_3R_Vz = xz_V(cs_3R_idx);
cs_3R = [cs_3R_P cs_3R_T cs_3R_M cs_3R_Vy cs_3R_Vz];
% Cross Section 4
[~, cs_4L_idx] = closest(xy_x, L_E4);
cs_4L_idx = cs_4L_idx - idx;
cs_4L_P = xy_P(cs_4L_idx);
cs_4L_T = xy_T(cs_4L_idx)*1e3; % [Nmm]
cs_4L_M = M(cs_4L_idx)*1e3; % [Nmm]
cs_4L_Vy = xy_V(cs_4L_idx);
cs_4L_Vz = xz_V(cs_4L_idx);
cs_4L = [cs_4L_P cs_4L_T cs_4L_M cs_4L_Vy cs_4L_Vz];
[~, cs_4R_idx] = closest(xy_x, L_E4);
cs_4R_idx = cs_4R_idx + idx;
cs_4R_P = xy_P(cs_4R_idx);
cs_4R_T = xy_T(cs_4R_idx)*1e3; % [Nmm]
cs_4R_M = M(cs_4R_idx)*1e3; % [Nmm]
cs_4R_Vy = xy_V(cs_4R_idx);
cs_4R_Vz = xz_V(cs_4R_idx);
cs_4R = [cs_4R_P cs_4R_T cs_4R_M cs_4R_Vy cs_4R_Vz];
% Cross Section 5
[~, cs_5L_idx] = closest(xy_x, L_E5);
cs_5L_idx = cs_5L_idx - idx;
cs_5L_P = xy_P(cs_5L_idx);
cs_5L_T = xy_T(cs_5L_idx)*1e3; % [Nmm]
cs_5L_M = M(cs_5L_idx)*1e3; % [Nmm]
cs_5L_Vy = xy_V(cs_5L_idx);
cs_5L_Vz = xz_V(cs_5L_idx);
cs_5L = [cs_5L_P cs_5L_T cs_5L_M cs_5L_Vy cs_5L_Vz];
[~, cs_5R_idx] = closest(xy_x, L_E5);
cs_5R_idx = cs_5R_idx + idx;
cs_5R_P = xy_P(cs_5R_idx);
cs_5R_T = xy_T(cs_5R_idx)*1e3; % [Nmm]
cs_5R_M = M(cs_5R_idx)*1e3; % [Nmm]
cs_5R_Vy = xy_V(cs_5R_idx);
cs_5R_Vz = xz_V(cs_5R_idx);
cs_5R = [cs_5R_P cs_5R_T cs_5R_M cs_5R_Vy cs_5R_Vz];
% Cross Section 6
[~, cs_6L_idx] = closest(xy_x, L_E6);
cs_6L_idx = cs_6L_idx - idx;
cs_6L_P = xy_P(cs_6L_idx);
cs_6L_T = xy_T(cs_6L_idx)*1e3; % [Nmm]
cs_6L_M = M(cs_6L_idx)*1e3; % [Nmm]
cs_6L_Vy = xy_V(cs_6L_idx);
cs_6L_Vz = xz_V(cs_6L_idx);
cs_6L = [cs_6L_P cs_6L_T cs_6L_M cs_6L_Vy cs_6L_Vz];
[~, cs_6R_idx] = closest(xy_x, L_E6);
cs_6R_idx = cs_6R_idx + idx;
cs_6R_P = xy_P(cs_6R_idx);
cs_6R_T = xy_T(cs_6R_idx)*1e3; % [Nmm]
cs_6R_M = M(cs_6R_idx)*1e3; % [Nmm]
cs_6R_Vy = xy_V(cs_6R_idx);
cs_6R_Vz = xz_V(cs_6R_idx);
cs_6R = [cs_6R_P cs_6R_T cs_6R_M cs_6R_Vy cs_6R_Vz];
%%%% For Fatigue %%%%
%%%% For Bearings %%%%
[~, E_idx] = closest(xy_x, 0);
E_Fa = xy_P(E_idx); % [N] Axial Force
E_Fr = sqrt( xy_V(E_idx)^2 + xz_V(E_idx)^2 ); % [N] Radial Force
[~, D_idx] = closest(xy_x, L_ED);
D_Fa = xy_P(D_idx); % [N] Axial Force
D_Fr = sqrt( xy_V(D_idx)^2 + xz_V(D_idx)^2 ); % [N] Radial Force
%%%% For Bearings %%%%

% Export data w/o figures (https://stackoverflow.com/questions/
                            % 45560181/avoid-saving-of-graphics-in-matlab)
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save(fullfile(export_import, "loadingDiagram_shaft2.mat"), saveVars{:})

%% Length sanity check
lW = 3;

if (~disablePlotting) && (~disablePlot2)
    figHandle = 12;
    figure(figHandle)
    hold on
    plot( [0 L_ED],  [0  0], 'k', 'LineWidth', lW)
    plot( [0 L_EG2], [2  2], 'LineWidth', lW)
    plot( [0 L_EG3], [-1 -1], 'LineWidth', lW)
    plot( [L_EG3 (L_EG3 + L_G3G2)], [-1 -1], 'LineWidth', lW)
    plot( [L_EG2 (L_EG2 + L_G2D)], [-1 -1], 'LineWidth', lW)
    plot( [L_EG3 (L_EG3 + L_G3D)], [-2 -2], 'LineWidth', lW)
    xlabel('Length [m]')
    legend('$L_{ED}$', '$L_{EG2}$', '$L_{EG3}$', '$L_{G3G2}$', ...
            '$L_{G2D}$', '$L_{G3D}$', ...
            'Location','southoutside', 'interpreter', 'latex')
    xlim( [ (-L_ED*0.1), (L_ED + L_ED*0.1) ] )
    ylim( [-5 5] )
    title('Shaft 2: One Directional Length', 'interpreter', 'latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAS413 Project: Loading Diagrams - Shaft 3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clearvars -except export_import;
load(fullfile(export_import, "loadingDiagram_common.mat"))

% Lengths
L_FG = b_F/2 + L_78 + b_s2 + L_45 + b_s1 + L_12 + b_G/2; % [m]
L_FG4 = b_F/2 + L_78 + b_s2/2; % [m]
L_G4G = L_FG - L_FG4; % [m]
L_FH = L_FG + L_GH; % [m]
L_F7 = b_F/2; % [m]
L_F8 = b_F/2 + L_78; % [m]
L_F9 = L_FG - b_G/2; % [m]

% Reaction forces @ bearings
F_Fz = (F_a4*r_G4 + F_r4*L_G4G) / L_FG; % [N]
F_Fy = (F_t4*L_G4G) / L_FG; % [N]
F_Fx = F_a4; % [N]
F_Gz = F_r4 - F_Fz; % [N]
F_Gy = F_t4 - F_Fy; % [N]

F_FG = (F_G4 *L_G4G )/(L_FG); %[N]
F_GG = F_G4 - F_FG; %[N]


%% XY - Plane

% Figure setup
figHandle = 7;
xPos = 10;
yPos = 3;

% Initialize XY Plots
xy_x = [];
xy_P = [];
xy_V = [];
xy_M = [];
xy_T = [];

% 0 < x < L_FG4
x = linspace(0, L_FG4, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Fx)]; % [N]
xy_V = [xy_V, ones(size(x)) * F_Fy]; % [N]
xy_M = [xy_M, (F_Fy*x)]; % [Nm]
xy_T = [xy_T, zeros(size(x) )]; % [Nm]

% L_FG4 < x < L_FG
x = linspace(L_FG4, L_FG, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x))* (F_a4 - F_Fx)]; % [N]
xy_V = [xy_V, ones(size(x)) * (F_Fy - F_t4)]; % [N]
xy_M = [xy_M, ( F_Fy*(x) - F_t4*(x - L_FG4) )]; % [Nm]
xy_T = [xy_T, ones(size(x))*(F_t4*r_G4)]; % [Nm]

% L_FG < x < L_FH
x = linspace(L_FG, L_FH, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Fx+F_a4)]; % [N]
xy_V = [xy_V, ones(size(x)) * (F_Fy - F_t4 + F_Gy)]; % [N]
xy_M = [xy_M, ( F_Fy*(x) - F_t4*(x - L_FG4) + F_Gy*(x - L_FG) )]; % [Nm]
xy_T = [xy_T, ones(size(x))*(F_t4*r_G4)]; % [Nm]

if (~disablePlotting) && (~disablePlot3)
    XYplaneFig = figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 3: XY - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_y(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_z(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')
    
    subplot(2,2,1)
    plotLD(xy_x,xy_P, colFill)
    subplot(2,2,2)
    plotLD(xy_x,xy_V, colFill)
    subplot(2,2,3)
    plotLD(xy_x,xy_M, colFill)
    subplot(2,2,4)
    plotLD(xy_x,xy_T, colFill)
end

%% XZ - Plane

% Figure setup
figHandle = 8;
xPos = 10;
yPos = 3;

% Initialize XZ Plots
xz_x = [];
xz_P = [];
xz_V = [];
xz_M = [];
xz_T = [];
xz_Mg = [];

% 0 < x < L_FG4
x = linspace(0, L_FG4, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) )]; % [Nm]
xz_T = [xz_T, zeros(size(x))]; % [Nm]
xz_Mg = [xz_Mg, F_FG*x];

% L_FG4 < x < L_FG
x = linspace(L_FG4, L_FG, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx + F_a4) ]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz - F_r4)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) - F_r4*(x - L_FG4) - F_a4*(r_G4) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * F_t4 *r_G4]; % [Nm]
xz_Mg = [xz_Mg, F_FG*x - F_G4*(x-L_FG4)];

% L_FG < x < L_FH
x = linspace(L_FG, L_FH, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx + F_a4)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz - F_r4 + F_Gz)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) - F_r4*(x - L_FG4) - F_a4*(r_G4) + ...
                F_Gz*(x - L_FG) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (F_t4*r_G4)]; % [Nm]
xz_Mg = [xz_Mg, F_FG*x - F_G4*(x-L_FG4)+F_GG*(x-L_FG)];

if (~disablePlotting) && (~disablePlot3)
    XZplaneFig = figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos wPlot hPlot]);
    sgtitle('\textbf{Shaft 3: XZ - Plane}', 'interpreter', 'latex')
    subplot(2,2,1)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Axial Force $P(x)$', 'Interpreter','latex')
    subplot(2,2,2)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[N]', 'interpreter', 'latex')
    title('Shear Force $V_z(x)$', 'Interpreter','latex')
    subplot(2,2,3)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Bending Moment $M_y(x)$', 'Interpreter','latex')
    subplot(2,2,4)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    title('Torque $T(x)$', 'Interpreter','latex')
    
    subplot(2,2,1)
    plotLD(xz_x,xz_P,colFill)
    subplot(2,2,2)
    plotLD(xz_x,xz_V,colFill)
    subplot(2,2,3)
    plotLD(xz_x,xz_M,colFill)
    subplot(2,2,4)
    plotLD(xz_x,xz_T,colFill)
end

%% Loading Diagram: Combined Bending Moment

M = sqrt(xz_M.^2 + xy_M.^2); % [Nm]
[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);

% Figure setup
figHandle = 9;
wPlotM = wPlot;
hPlotM = hPlot/2;

if (~disablePlotting) && (~disablePlot3)
    figure(figHandle);
    set(figHandle,'Units','Centimeter')
    set(figHandle,'Position',[xPos yPos+hPlotM/2 wPlotM hPlotM]);
    plotLD(xy_x, M, colFill)
    title('\textbf{Shaft 3: Combined Bending Moment}', ...
        'interpreter', 'latex', 'FontSize',fSize)
    xlabel('[m]', 'interpreter', 'latex')
    ylabel('[Nm]', 'interpreter', 'latex')
    xlim([0 L_FH])
    
    dashLineV(L, figHandle, 2, 2)
end

%% Export
idx = 2;
%%%% For Fatigue %%%%
% Cross Section 7
[~, cs_7L_idx] = closest(xy_x, L_F7);
cs_7L_idx = cs_7L_idx - idx;
cs_7L_P = xy_P(cs_7L_idx);
cs_7L_T = xy_T(cs_7L_idx)*1e3; % [Nmm]
cs_7L_M = M(cs_7L_idx)*1e3; % [Nmm]
cs_7L_Vy = xy_V(cs_7L_idx);
cs_7L_Vz = xz_V(cs_7L_idx);
cs_7L = [cs_7L_P cs_7L_T cs_7L_M cs_7L_Vy cs_7L_Vz];
[~, cs_7R_idx] = closest(xy_x, L_F7);
cs_7R_idx = cs_7R_idx + idx;
cs_7R_P = xy_P(cs_7R_idx);
cs_7R_T = xy_T(cs_7R_idx)*1e3; % [Nmm]
cs_7R_M = M(cs_7R_idx)*1e3; % [Nmm]
cs_7R_Vy = xy_V(cs_7R_idx);
cs_7R_Vz = xz_V(cs_7R_idx);
cs_7R = [cs_7R_P cs_7R_T cs_7R_M cs_7R_Vy cs_7R_Vz];
% Cross Section 8
[~, cs_8L_idx] = closest(xy_x, L_F8);
cs_8L_idx = cs_8L_idx - idx;
cs_8L_P = xy_P(cs_8L_idx);
cs_8L_T = xy_T(cs_8L_idx)*1e3; % [Nmm]
cs_8L_M = M(cs_8L_idx)*1e3; % [Nmm]
cs_8L_Vy = xy_V(cs_8L_idx);
cs_8L_Vz = xz_V(cs_8L_idx);
cs_8L = [cs_8L_P cs_8L_T cs_8L_M cs_8L_Vy cs_8L_Vz];
[~, cs_8R_idx] = closest(xy_x, L_F8);
cs_8R_idx = cs_8R_idx + idx;
cs_8R_P = xy_P(cs_8R_idx);
cs_8R_T = xy_T(cs_8R_idx)*1e3; % [Nmm]
cs_8R_M = M(cs_8R_idx)*1e3; % [Nmm]
cs_8R_Vy = xy_V(cs_8R_idx);
cs_8R_Vz = xz_V(cs_8R_idx);
cs_8R = [cs_8R_P cs_8R_T cs_8R_M cs_8R_Vy cs_8R_Vz];
% Cross Section 9
[~, cs_9L_idx] = closest(xy_x, L_F9);
cs_9L_idx = cs_9L_idx - idx;
cs_9L_P = xy_P(cs_9L_idx);
cs_9L_T = xy_T(cs_9L_idx)*1e3; % [Nmm]
cs_9L_M = M(cs_9L_idx)*1e3; % [Nmm]
cs_9L_Vy = xy_V(cs_9L_idx);
cs_9L_Vz = xz_V(cs_9L_idx);
cs_9L = [cs_9L_P cs_9L_T cs_9L_M cs_9L_Vy cs_9L_Vz];
[~, cs_9R_idx] = closest(xy_x, L_F9);
cs_9R_idx = cs_9R_idx + idx;
cs_9R_P = xy_P(cs_9R_idx);
cs_9R_T = xy_T(cs_9R_idx)*1e3; % [Nmm]
cs_9R_M = M(cs_9R_idx)*1e3; % [Nmm]
cs_9R_Vy = xy_V(cs_9R_idx);
cs_9R_Vz = xz_V(cs_9R_idx);
cs_9R = [cs_9R_P cs_9R_T cs_9R_M cs_9R_Vy cs_9R_Vz];
% Cross Section H
[~, cs_H_idx] = closest(xy_x, L_FH);
cs_H_P = xy_P(cs_H_idx);
cs_H_T = xy_T(cs_H_idx)*1e3; % [Nmm]
cs_H_M = M(cs_H_idx)*1e3; % [Nmm]
cs_H_Vy = xy_V(cs_H_idx);
cs_H_Vz = xz_V(cs_H_idx);
cs_H = [cs_H_P cs_H_T cs_H_M cs_H_Vy cs_H_Vz];
%%%% For Fatigue %%%%
%%%% For Bearings %%%%
[~, F_idx] = closest(xy_x, 0);
F_Fa = xy_P(F_idx); % [N] Axial Force
F_Fr = sqrt( xy_V(F_idx)^2 + xz_V(F_idx)^2 ); % [N] Radial Force
[~, G_idx] = closest(xy_x, L_FG);
G_Fa = xy_P(G_idx); % [N] Axial Force
G_Fr = sqrt( xy_V(G_idx)^2 + xz_V(G_idx)^2 ); % [N] Radial Force
%%%% For Bearings %%%%

% Export data w/o figures (https://stackoverflow.com/questions/
                            % 45560181/avoid-saving-of-graphics-in-matlab)
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save(fullfile(export_import, "loadingDiagram_shaft3.mat"), saveVars{:})

%% Length sanity check
lW = 3;

if (~disablePlotting) && (~disablePlot3)
    figHandle = 10;
    figure(figHandle)
    hold on
    plot( [0 L_FH],  [0  0], 'k', 'LineWidth', lW)
    plot( [0 L_FG], [1  1], 'LineWidth', lW)
    plot( [0 L_FG4], [-1 -1], 'LineWidth', lW)
    plot( [L_FG4 (L_FG4 + L_G4G)], [-1 -1], 'LineWidth', lW)
    plot( [L_FG (L_FG + L_GH)], [-1 -1], 'LineWidth', lW)
    xlabel('Length [m]')
    legend('$L_{FH}$', '$L_{FG}$', '$L_{FG4}$', '$L_{G4G}$', '$L_{GH}$', ...
            'Location','southoutside', 'interpreter', 'latex')
    xlim( [ (-L_FH*0.1), (L_FH + L_FH*0.1) ] )
    ylim( [-5 5] )
    title('Shaft 3: One Directional Length', 'interpreter', 'latex')
end



%% Common functions for loading diagrams

function plotLD(x, f, colFill)
    % Plot f(x) as a black line with shaded area 
    % between f(x) and abscissa

    % Aesthetics :sparkly_glitter:
    colLine = 'k';
    lineW = 2;

    % Filled Area
    X = [x, fliplr(x)];
    Y = [f, zeros(size(f))];

    % Actual Plotting
    hold on
    fill(X, Y, colFill, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(x, f,'Color', colLine','LineWidth', lineW)
    yline(0, 'k')

    % Detect Jumps and Draw Line
    jumpTol = 1;
    jumps = find( abs( diff(f) ) > jumpTol );
    for i = 1 : length(jumps)
        xJump = x(jumps(i) + 1);
        fPreJump = f(jumps(i));
        fPostJump = f(jumps(i));
        plot( [ xJump xJump ], [ fPreJump, fPostJump] , ...
              colLine, 'LineWidth', lineW)
    end

end

function dashLineV(L, figNr, spRow, spCol)
    % Add dashed line to all figures
    for i = (figNr-2) : figNr
        if i < figNr
            figure(i)
            for row = 1 : spRow
                for col = 1 : spCol
                    spIdx = (row-1) * spCol + col;
                    subplot(spRow,spCol,spIdx)
                    xline(L, 'm--', 'LineWidth', 1)
                end
            end
        else
            figure(i)
            xline(L, 'm--', 'LineWidth', 1)
        end
    end
end