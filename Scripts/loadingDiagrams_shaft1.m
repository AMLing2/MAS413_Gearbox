clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 1 %
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
L_AB  = 0.05; % [m]
    % Bearing widths
b_B = 30e-3; % [m] catalogue circa 16 - 47 [mm] <-- WIP
b_C = b_B; % [m]

% Import from Gear Sizing
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', 'b_s1', 'b_s2')
    % Convert from Gear Sizing
r_G1 = d_g1/2 * 1e-3; % [m]
r_G2 = d_g2/2 * 1e-3; % [m]
r_G3 = d_g3/2 * 1e-3; % [m]
r_G4 = d_g4/2 * 1e-3; % [m]
b_s1 = b_s1 * 1e-3; % [m]
b_s2 = b_s2 * 1e-3; % [m]

% Calculated values
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
T_M = P_1 / omega_1; % [Nm]
    % Gear Forces
F_t1 = T_M / r_G1; % [N]
F_a1 = F_t1 * tand(beta); % [N]
F_r1 = F_t1 * tand(alpha)/cosd(beta); % [N]
    % Lenghts
L_BC  = b_B/2 + L_78 + b_s2 + L_45 + b_s1 + L_12; % [m]
L_BG1 = b_B/2 + L_78 + b_s2 + L_45 + b_s1/2; % [m]
L_G1C = L_BC - L_BG1; % [m]
L_AG1 = L_AB + L_BG1; % [m]
L_AC = L_AG1 + L_G1C; % [m]
    % For Reaction forces @ bearings
F_By = F_t1*L_G1C/L_BC; % [N]
F_Bz = (F_r1*L_G1C - F_a1*r_G1)/L_BC; % [N]


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
xy_x = [xy_x, x];
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, zeros(size(x))]; % [N]
xy_M = [xy_M, zeros(size(x))]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M ]; % [Nm]

% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, zeros(size(x))]; % [N]
xy_V = [xy_V, ones(size(x)) * F_By]; % [N]
xy_M = [xy_M, ( F_By*(x - L_AB) )]; % [Nm]
xy_T = [xy_T, ones(size(x)) * T_M]; % [Nm]

% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xy_x = [xy_x, x];
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
xz_M = [xz_M, ( - F_Bz*(x - L_AB) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]


% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_a1)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r1 - F_Bz)]; % [N]
xz_M = [xz_M, ( F_r1*(x - L_AG1) - F_Bz*(x - L_AB) - F_a1*(r_G1) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (T_M - F_t1*r_G1)]; % [Nm]


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

%% Shaft deflection calculations

%E-modul
E =210; %Emodul - Change 

I_shaft =[];

%Diameters of shaft
d_c  =1; %[mm]
d_12 =1; %[mm]
d_S1 =1; %[mm]

%Lengths of shaft
b_G1 = 1; %[mm]

%Calculate I for the different intervals

% 0 < x < L_AB
x = linspace(0, L_AB, resolution); 
I_shaft = [I_shaft, ones(size(x))]; % Moment of area
% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
I_shaft = [I_shaft, ones(size(x))]; % Moment of area
% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
I_shaft = [I_shaft, ones(size(x))];



function [Iy, Iz, Ix] = secondMomentAreaCyl(D)
    % bruk: [Iy, Iz, Ix] = secondMomentAreaCyl(D)
    % D - diameter til sylinderen
    % Iy, Iz - Arealmoment om y- og z-aksene
    % Ix - Arealmoment om x-aksen (langs sylinderens lengde)
    R=D/2;
    Iy = (pi / 4) * R^4;  % Om y-aksen
    Iz = Iy;              % Om z-aksen 
    Ix = (pi / 2) * R^4;  % Om x-aksen 
end
