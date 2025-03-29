clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

% Masses of gears
g = 9.81; %[m/s^2]
m_G2 = 0.52675; %[kg]
m_G3 = 1; %[kg]

% Diameters of shaft
d_D = 0.01; % [m]
d_E = 0.015; % [m]
d_S21 = 0.02; % [m]
d_S22 = 0.02; %[m]
d_45 = 0.02; %[m]

% Lengths of gears
b_s1 = 0.005; % [m]
b_s2 = 0.005; % [m]

%Lengths of bearings
L_BE = 0.005; %[m]
L_BD = 0.005; %[m]

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

% Chosen Parameters
L_12 = 5e-3; % [m]
L_45 = 5e-3; % [m]
L_78 = 5e-3; % [m]
    % Bearing widths
b_D = 30e-3; % [m] catalogue circa 16 - 47 [mm] <-- WIP
b_E = b_D; % [m]
% eta = 0.96; % [-] Stage efficiency "finely worked teeth & good lubrication"
eta = 1.00; % [-] Ideal Stages
 
% Import from Gear Sizing
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', 'b_s1', 'b_s2', 'i_tot', 'i_s1', 'E_mat')
    % Convert from Gear Sizing
r_G1 = d_g1/2 * 1e-3; % [m]
r_G2 = d_g2/2 * 1e-3; % [m]
r_G3 = d_g3/2 * 1e-3; % [m]
r_G4 = d_g4/2 * 1e-3; % [m]
b_s1 = b_s1 * 1e-3; % [m]
b_s2 = b_s2 * 1e-3; % [m]
E = E_mat*1e6;%[Pa]

% Calculated Values
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
n_out = (n_1/i_tot); % [RPM]
omega_out = n_out * 2*pi / 60; % [rad/sec]
eta_tot = eta^2; % [-] Squared because there are two stages
P_out = P_1*eta_tot; % [W]
T_M   = P_1/omega_1; % [Nm]
T_out = P_out/omega_out; %[Nm] 
    % Lengths
L_ED = b_E/2 + L_78 + b_s2 + L_45 + b_s1 + L_12 + b_D/2; % [m]
L_EG3 = b_E/2 + L_78 + b_s2/2; % [m]
L_EG2 = L_EG3 + b_s2/2 + L_45 + b_s1/2; % [m]
L_G3G2 = L_EG2 - L_EG3; % [m]
L_G3D = L_ED - L_EG3; % [m]
L_G2D = L_ED - L_EG2; % [m]
L_E3 = L_ED - b_E/2; % [m]
L_E4 = L_EG2 - b_s1/2; % [m]
L_E5 = L_EG3 + b_s2/2; % [m]
L_E6 = b_E/2; % [m]
    % Gear 2 forces
F_t2 = T_M/r_G1; % [N]
F_a2 = F_t2 * tand(beta); % [N]
F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]
F_G2 = (m_G2*g);  %[N]

    % Gear 3 forces
F_t3 = T_out / r_G4; % [N]
F_a3 = F_t3 * tand(beta); % [N]
F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]
F_G3 = (m_G3*g); %[N]

    % Reaction forces @ bearings
F_Ex = F_a2 - F_a3;
F_Ez = (F_r3*L_G3D - F_r2*L_G2D - F_a3*r_G3 - F_a2*r_G2)/L_ED;
F_Ey = (F_t3*L_G3D + F_t2*L_G2D)/L_ED;

F_EG = ( (F_G3*L_G3D) + ( F_G2*L_G2D) ) / (L_ED);

 
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
title('Torque $T(x)$', 'Interpreter','latex')

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

%% Export

%%%% For Fatigue %%%%
% Cross Section 3
[~, cs_3_idx] = closest(xy_x, L_E3);
cs_3_P = xy_P(cs_3_idx);
cs_3_T = xy_T(cs_3_idx)*1e3; % [Nmm]
cs_3_M = M(cs_3_idx)*1e3; % [Nmm]
cs_3_Vy = xy_V(cs_3_idx);
cs_3_Vz = xz_V(cs_3_idx);
cs_3 = [cs_3_P cs_3_T cs_3_M cs_3_Vy cs_3_Vz];
% Cross Section 4
[~, cs_4_idx] = closest(xy_x, L_E4);
cs_4_P = xy_P(cs_4_idx);
cs_4_T = xy_T(cs_4_idx)*1e3; % [Nmm]
cs_4_M = M(cs_4_idx)*1e3; % [Nmm]
cs_4_Vy = xy_V(cs_4_idx);
cs_4_Vz = xz_V(cs_4_idx);
cs_4 = [cs_4_P cs_4_T cs_4_M cs_4_Vy cs_4_Vz];
% Cross Section 5
[~, cs_5_idx] = closest(xy_x, L_E5);
cs_5_P = xy_P(cs_5_idx);
cs_5_T = xy_T(cs_5_idx)*1e3; % [Nmm]
cs_5_M = M(cs_5_idx)*1e3; % [Nmm]
cs_5_Vy = xy_V(cs_5_idx);
cs_5_Vz = xz_V(cs_5_idx);
cs_5 = [cs_5_P cs_5_T cs_5_M cs_5_Vy cs_5_Vz];
% Cross Section 6
[~, cs_6_idx] = closest(xy_x, L_E6);
cs_6_P = xy_P(cs_6_idx);
cs_6_T = xy_T(cs_6_idx)*1e3; % [Nmm]
cs_6_M = M(cs_6_idx)*1e3; % [Nmm]
cs_6_Vy = xy_V(cs_6_idx);
cs_6_Vz = xz_V(cs_6_idx);
cs_6 = [cs_6_P cs_6_T cs_6_M cs_6_Vy cs_6_Vz];
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
save("loadingDiagram_shaft2.mat", saveVars{:})

%% Length sanity check
lW = 3;

figure(99)
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
title('One Directional Length', 'interpreter', 'latex')