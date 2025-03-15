clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 3 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

% Common Plotting Constants
colFill = [0.7765 0.9176 0.9843];
resolution = 100;
wPlot = 22;
hPlot = 16;
fSize = 16;

% Given information
n_1 = 1450; % [RPM]
P_1 = 12.5e3; % [W]
i_tot = 17.3;
alpha = 20; % [degrees] Helix Angle
beta = 15;  % [degrees] Pressure Angle

% Chosen Parameters
L_12 = 5e-3; % [m]
L_45 = 5e-3; % [m]
L_78 = 5e-3; % [m]
L_GH = 0.05; % [m]
    % Bearing widths
b_F = 30e-3; % [m] catalogue circa 16 - 47 [mm] <-- WIP
b_G = b_F; % [m]
% eta = 0.96; % [-] Stage efficiency "finely worked teeth & good lubrication"
eta = 1.00; % [-] Ideal Stages

% Import from Gear Sizing
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', 'b_s1', 'b_s2')
    % Convert from Gear Sizing
r_G1 = d_g1/2 * 1e-3; % [m]
r_G2 = d_g2/2 * 1e-3; % [m]
r_G3 = d_g3/2 * 1e-3; % [m]
r_G4 = d_g4/2 * 1e-3; % [m]
b_s1 = b_s1 * 1e-3; % [m]
b_s2 = b_s2 * 1e-3; % [m]

% Calculated Values
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
n_out = (n_1/i_tot); % [RPM]
omega_out = n_out * 2*pi / 60; % [rad/sec]
eta_tot = eta^2; % (Squared efficiency because there are two stages)
P_out = P_1*eta_tot; % [W]
T_M = P_1/omega_1; % [Nm]
T_out = P_out/omega_out; % [Nm] 
    % Lengths
L_FG = b_F/2 + L_78 + b_s2 + L_45 + b_s1 + L_12 + b_G/2; % [m]
L_FG4 = b_F/2 + L_78 + b_s2/2; % [m]
L_G4G = L_FG - L_FG4; % [m]
L_FH = L_FG + L_GH; % [m]
    % Gear 3 forces
F_t4 = (T_out*i_tot) / r_G4; % [N]
F_a4 = F_t4 * tand(beta); % [N]
F_r4 = F_t4 * tand(alpha)/cosd(beta); % [N]
    % Reaction forces @ bearings
F_Fz = (F_a4*r_G4 + F_r4*L_G4G) / L_FG; % [N]
F_Fy = (F_t4*L_G4G) / L_FG; % [N]
F_Fx = F_a4; % [N]
F_Gz = F_r4 - F_Fz; % [N]
F_Gy = F_t4 - F_Fy; % [N]


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

% 0 < x < L_FG4
x = linspace(0, L_FG4, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) )]; % [Nm]
xz_T = [xz_T, zeros(size(x))]; % [Nm]

% L_FG4 < x < L_FG
x = linspace(L_FG4, L_FG, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx + F_a4) ]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz - F_r4)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) - F_r4*(x - L_FG4) - F_a4*(r_G4) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * F_t4 *r_G4]; % [Nm]

% L_FG < x < L_FH
x = linspace(L_FG, L_FH, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_Fx + F_a4)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_Fz - F_r4 + F_Gz)]; % [N]
xz_M = [xz_M, ( F_Fz*(x) - F_r4*(x - L_FG4) - F_a4*(r_G4) + ...
                F_Gz*(x - L_FG) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (F_t4*r_G4)]; % [Nm]

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
xlim([0 L_FH])

[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);
dashLineV(L, 3, 2, 2)

%% Length sanity check
lW = 3;

figure(99)
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
title('One Directional Length', 'interpreter', 'latex')