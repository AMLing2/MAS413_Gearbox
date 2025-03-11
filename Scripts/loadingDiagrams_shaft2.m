clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

% Common Plotting Constants
colFill = [0.7765 0.9176 0.9843];
resolution = 100;
wPlot = 22;
hPlot = 16;
fSize = 16;

% Given information from project description - Needs adjustment of values
i_tot = 17.3; % [-]
i_1 = 1; % [-]
i_2 = 2; % [-]


n_in = 1450; % [RPM]
omega_in = n_in * 2*pi / 60; % [rad/sec]
n_out = (n_in/i_tot); %[RPM]
omega_out = n_out * 2*pi / 60; %[rad/sec]
 
eta = 0.96; %[-] Efficiency
eta_tot = eta^2;

P_in = 12.5e3; % [W]
P_out = P_in*eta_tot; %[W] (Squared efficiency because there are two stages)
T_M = P_in/omega_in; %[Nm]
T_out = P_out/omega_out; %[Nm] 
 
alpha = 20; % [degrees] Helix Angle
beta = 15;  % [degrees] Pressure Angle

% Lengths on shaft - Needs adjustment of value
L_EG3  = 0.05; % [m]
L_G3G2 = 0.10; % [m]
L_G2D = 0.15; % [m]
L_tot = L_EG3 + L_G3G2 + L_G2D; % [m]
L_G3D = L_G3G2+L_G2D; % [m]
 
% Radius of gears - Calculated Elsewhere - Needs adjustment of value
r_G1 = 0.25; % [m]
r_G2 = 0.25; % [m]
r_G3 = 0.25; % [m]
r_G4 = 0.25; % [m]

%% Calculations
% Calculated lengths & Torques
L_E_G2 = L_EG3 + L_G3G2; % [m]
 
%Gear 2 forces
F_t2 = T_M/r_G1; % [N]
F_a2 = F_t2 * tand(beta); % [N]
F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]

% Gear 3 forces
F_t3 = (T_out*i_tot) / r_G4; % [N]
F_a3 = F_t3 * tand(beta); % [N]
F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]

%% Reaction forces @ bearings
L_ED = L_tot;

F_Ex = F_a2-F_a3;

F_Ez = (F_r3*L_G3D - F_r2*L_G2D - F_a3*r_G3 - F_a2*r_G2)/L_ED;

F_Ey= (F_t3*L_G3D + F_t2*L_G2D)/L_ED;

 
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
title('Axial Torque $T(x)$', 'Interpreter','latex')

 
% 0 < x < L_E_G3
x = linspace(0, L_EG3, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Ex)]; % [N] %%% changed sign of F_Ex
xy_V = [xy_V, ones(size(x)) * (-F_Ey)]; % [N]
xy_M = [xy_M, (-F_Ey * x)]; % [Nm]
xy_T = [xy_T, zeros(size(x))]; % [Nm]


 
% L_E_G3 < x < L_G3_G2
x = linspace(L_EG3, L_G3G2, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x))* -(F_Ex+F_a3)]; % [N]
xy_V = [xy_V, ones(size(x)) * (-F_Ey+F_t3)]; % [N] %%% changed sign of F_t3
xy_M = [xy_M,-F_Ey*x + F_t3*(x - L_EG3)]; % [Nm]
xy_T = [xy_T,ones(size(x)) * F_t3*r_G3]; % [Nm]
 
 
% L_G3_G2 < x < L_G2_D 
x = linspace(L_G3G2, L_G2D, resolution);
xy_x = [xy_x, x];
xy_P = [xy_P, ones(size(x)) * (-F_Ex - F_a2 + F_a3)]; % [N] %%% changed sign of F_Ex and F_a3
xy_V = [xy_V,ones(size(x)) * (F_t2 - F_Ey + F_t3)]; % [N]
xy_M = [xy_M,F_t2 * (x- L_E_G2) - F_Ey*x + F_t3 * (x - L_EG3)]; % [Nm]
xy_T = [xy_T,ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2)]; % [Nm]
 
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
title('Axial Torque $T(x)$', 'Interpreter','latex')
 
 
 % 0 < x < L_E_G3
 
 x = linspace(0, L_EG3, resolution);
 xz_x = [xz_x, x];
 xz_P = [xz_P,ones(size(x))*(-F_Ex)]; % [N]
 xz_V = [xz_V,ones(size(x))*(-F_Ez)]; % [N]
 xz_M = [xz_M,ones(size(x)).* (x *-(F_Ez))]; % [Nm] %%% changed -(F_ex) to -(F_ez)
 xz_T = [xz_T,zeros(size(x))]; % [Nm]
 

% L_E_G3 < x < L_G3_G2
x = linspace(L_EG3, L_G3G2, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P,ones(size(x))* -(F_Ex+F_a3)]; % [N]
xz_V = [xz_V,ones(size(x)) * (F_r3-F_Ez)]; % [N]
xz_M = [xz_M,F_Ez*x - F_r3*(x - L_EG3)]; % [Nm] %%%%%%%%%%%%%%% missing F_a3*r_G3
xz_T = [xz_T,ones(size(x)) * F_t3*r_G3]; % [Nm]
 
% L_G3_G2 < x < L_G2_D

x = linspace(L_G3G2, L_G2D, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P,ones(size(x)) * (F_a2 - F_a3 - F_Ex)]; % [N]
xz_V = [xz_V,ones(size(x)) * (F_r3 - F_Ez - F_r2)]; % [N] %%%%%%%%%%%%%  wrong
xz_M = [xz_M,F_Ez * x - F_r3*(x - L_EG3) + F_a3 * r_G3 * (x - L_EG3) + F_a2*r_G2*(x-L_G3G2)]; % [Nm] %%%%% Nm^2 for F_a2...
xz_T = [xz_T, ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2)]; % [Nm]
 
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
xlim([0 L_ED])

[M_max, M_max_idx] = max(M);
L = xy_x(M_max_idx);
dashLineV(L, 3, 2, 2)