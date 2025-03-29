clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

% Chosen Parameters
L_12 = 5e-3; % [m]
L_45 = 5e-3; % [m]
L_78 = 5e-3; % [m]
L_AB  = 0.05; % [m]

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

% Import from Bearings
load('bearings.mat', 'b_B', 'b_C') % [mm]
b_B = b_B / 1000; % [m]
b_C = b_C / 1000; % [m]

% Diameters of shaft
d_C   = 0.010; % [m]
d_12  = 0.015; % [m]
d_B   = 0.011; % [m]
d_S1  = 0.013; % [m]

% Import from Gear Sizing
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', ...
    'b_s1', 'b_s2', 'i_tot', 'E_mat', 'mass_g1') % [mm], [MPa], [kg]
    
% Convert from Gear Sizing
r_G1 = d_g1/2 * 1e-3; % [m]
r_G2 = d_g2/2 * 1e-3; % [m]
r_G3 = d_g3/2 * 1e-3; % [m]
r_G4 = d_g4/2 * 1e-3; % [m]
b_s1 = b_s1 * 1e-3; % [m]
b_s2 = b_s2 * 1e-3; % [m]
E = E_mat*1e6; % [Pa]

% Calculated values
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
T_M = P_1 / omega_1; % [Nm]
F_G1 = mass_g1*g;
% Gear Forces
F_t1 = T_M / r_G1; % [N]
F_a1 = F_t1 * tand(beta); % [N]
F_r1 = F_t1 * tand(alpha)/cosd(beta); % [N]
% Lengths
L_BC  = b_B/2 + L_78 + b_s2 + L_45 + b_s1 + L_12; % [m]
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

XYplaneFig = figure(figHandle);
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

subplot(2,2,1)
plotLD(xy_x,xy_P, colFill)
subplot(2,2,2)
plotLD(xy_x,xy_V, colFill)
subplot(2,2,3)
plotLD(xy_x,xy_M, colFill)
subplot(2,2,4)
plotLD(xy_x,xy_T, colFill)

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

XZplaneFig = figure(figHandle);
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

subplot(2,2,1)
plotLD(xz_x,xz_P,colFill)
subplot(2,2,2)
plotLD(xz_x,xz_V,colFill)
subplot(2,2,3)
plotLD(xz_x,xz_M,colFill)
subplot(2,2,4)
plotLD(xz_x,xz_T,colFill)


%% Loading Diagram: Combined Bending Moment

M = sqrt(xz_M.^2 + xy_M.^2); % [Nm]

% Figure setup
figHandle = 3;
wPlotM = wPlot;
hPlotM = hPlot/2;

comboMomentFig = figure(figHandle);
set(figHandle,'Units','Centimeter')
set(figHandle,'Position',[xPos yPos+hPlotM/2 wPlotM hPlotM]);
plotLD(xy_x, M, colFill)
title('Combined Moment', 'interpreter', 'latex', 'FontSize',fSize)
xlabel('[m]', 'interpreter', 'latex')
ylabel('[Nm]', 'interpreter', 'latex')
xlim([0 L_AC])

[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);
dashLineV(L, 3, 2, 2)


%% Export

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
[~, cs_0_idx] = closest(xy_x, L_A0);
cs_0_P = xy_P(cs_0_idx); % [N]
cs_0_T = xy_T(cs_0_idx)*1e3; % [Nmm]
cs_0_M = M(cs_0_idx)*1e3; % [Nmm]
cs_0_Vy = xy_V(cs_0_idx); % [N]
cs_0_Vz = xz_V(cs_0_idx); % [N]
cs_0 = [cs_0_P cs_0_T cs_0_M cs_0_Vy cs_0_Vz];
% Cross Section 1
[~, cs_1_idx] = closest(xy_x, L_A1);
cs_1_P = xy_P(cs_1_idx); % [N]
cs_1_T = xy_T(cs_1_idx)*1e3; % [Nmm]
cs_1_M = M(cs_1_idx)*1e3; % [Nmm]
cs_1_Vy = xy_V(cs_1_idx); % [N]
cs_1_Vz = xz_V(cs_1_idx); % [N]
cs_1 = [cs_1_P cs_1_T cs_1_M cs_1_Vy cs_1_Vz];
% Cross Section 2
[~, cs_2_idx] = closest(xy_x, L_A2);
cs_2_P = xy_P(cs_2_idx); % [N]
cs_2_T = xy_T(cs_2_idx)*1e3; % [Nmm]
cs_2_M = M(cs_2_idx)*1e3; % [Nmm]
cs_2_Vy = xy_V(cs_2_idx); % [N]
cs_2_Vz = xz_V(cs_2_idx); % [N]
cs_2 = [cs_2_P cs_2_T cs_2_M cs_2_Vy cs_2_Vz];
%%%% For Fatigue %%%%
%%%% For Bearings %%%%
[~, B_idx] = closest(xy_x, L_AB);
B_Fa = xy_P(B_idx); % [N] Axial Force
B_Fr = sqrt( xy_V(B_idx)^2 + xz_V(B_idx)^2 ); % [N] Radial Force
[~, C_idx] = closest(xy_x, L_AC);
C_Fa = xy_P(C_idx); % [N] Axial Force
C_Fr = sqrt( xy_V(C_idx)^2 + xz_V(C_idx)^2 ); % [N] Radial Force
%%%% For Bearings %%%%

% Export data w/o figures (https://stackoverflow.com/questions/
                            % 45560181/avoid-saving-of-graphics-in-matlab)
varData = whos;
saveIndex = cellfun(@isempty, regexp({varData.class}, 'matlab.(graphics|ui)'));
saveVars = {varData(saveIndex).name};
save("loadingDiagram_shaft1.mat", saveVars{:})

%% Length sanity check
lW = 3;

figure(99)
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
title('One Directional Length', 'interpreter', 'latex')