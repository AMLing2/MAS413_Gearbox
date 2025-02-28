clc; clear; close all;
% MAS413 Project: Mechanics of Materials - Shaft 1

%% Constants

% Common Plotting Constants
colFill = [0.7765 0.9176 0.9843];
resolution = 100;
wPlot = 22;
hPlot = 16;

% Given information
n_1 = 1450; % [RPM]
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
P_1 = 12.5e3; % [W]
i_tot = 17.3;
alpha = 20; % [degrees] Helix Angle
beta = 15; % [degrees] Pressure Angle

% Chosen Parameters
L_AB  = 0.05; % [m]
L_BG1 = 0.10; % [m]
L_G1C = 0.15; % [m]

% Calculated Elsewhere
r_G1 = 0.25; % [m]

% Calculated values
L_AG1 = L_AB + L_BG1; % [m]
L_AC = L_AG1 + L_G1C; % [m]
T_M = P_1 / omega_1; % [Nm]
F_t1 = T_M / r_G1; % [N]
F_a1 = F_t1 * tand(beta); % [N]
F_r1 = F_t1 * tand(alpha)/cosd(beta); % [N]
    
% For Reaction forces @ bearings
L_BC = L_BG1 + L_G1C; % [m]

F_By = F_t1*L_G1C/L_BC; % [N]
F_Bz = F_r1*L_G1C/L_BC; % [N]


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
title('Axial Torque $T(x)$', 'Interpreter','latex')

% 0 < x < L_AB
x = linspace(0, L_AB, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, zeros(size(x))]; % [N]
xy_M = [xy_M, zeros(size(x))]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M ]; % [Nm]

% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, ones(size(x)) * (-F_By)]; % [N]
xy_M = [xy_M, F_By * (x - L_AB)]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M]; % [Nm]

% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_a1)]; % [N]
xy_V = [xy_V, ones(size(x)) * (F_By - F_t1)]; % [N]
xy_M = [xy_M, F_By * (x - L_AB) - F_t1 * (x - L_AG1)]; % [Nm]
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
title('Axial Torque $T(x)$', 'Interpreter','latex')

% 0 < x < L_AB
x = linspace(0, L_AB, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, zeros(size(x))]; % [N]
xz_M = [xz_M, zeros(size(x))]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]

% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, ones(size(x)) * (-F_Bz)]; % [N]
xz_M = [xz_M, - F_Bz * (x - L_AB)]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]

% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_a1)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r1 - F_Bz)]; % [N]
xz_M = [xz_M, F_r1 * (x - L_AG1) - F_Bz * (x - L_AB) - F_a1 * r_G1]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (T_M - F_t1*r_G1)]; % [Nm]

subplot(2,2,1)
plotLD(xz_x,xz_P,colFill)
subplot(2,2,2)
plotLD(xz_x,xz_V,colFill)
subplot(2,2,3)
plotLD(xz_x,xz_M,colFill)
subplot(2,2,4)
plotLD(xz_x,xz_T,colFill)