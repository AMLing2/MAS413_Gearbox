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
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', 'b_s1', 'b_s2', 'i_tot')
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
    % Gear 2 forces
F_t2 = T_M/r_G1; % [N]
F_a2 = F_t2 * tand(beta); % [N]
F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]
    % Gear 3 forces
F_t3 = T_out / r_G4; % [N]
F_a3 = F_t3 * tand(beta); % [N]
F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]
    % Reaction forces @ bearings
F_Ex = F_a2 - F_a3;
F_Ez = (F_r3*L_G3D - F_r2*L_G2D - F_a3*r_G3 - F_a2*r_G2)/L_ED;
F_Ey = (F_t3*L_G3D + F_t2*L_G2D)/L_ED;

 
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

% L_EG3 < x < L_EG2
x = linspace(L_EG3, L_EG2, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x))* ( - F_Ex - F_a3)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r3 - F_Ez)]; % [N]
xz_M = [xz_M, ( -F_Ez*(x) + F_r3*(x - L_EG3) - F_a3*(r_G3) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * F_t3*r_G3]; % [Nm]

% L_EG2 < x < L_ED
x = linspace(L_EG2, L_ED, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (F_a2 - F_a3 - F_Ex)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r3 - F_Ez - F_r2)]; % [N]
xz_M = [xz_M, ( -F_Ez*(x) + F_r3*(x - L_EG3) - F_r2*(x - L_EG2) - ...
                 F_a2*(r_G2) - F_a3*(r_G3) )]; % [Nm]
xz_T = [xz_T,  ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2)]; % [Nm]
 
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

save("loadingDiagram_shaft2.mat")

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


%% Shaft deflection calculations

% Visuals
lwDeflection = 2;
sizeDeflectionText = 16;
ok = 10;

res = 300;

%Initialize arrays:
I_shaft = [];
theta = zeros(1, res);
theta_corrected = [];
delta = zeros(1, res);
delta_corrected = [];

E = 210e9; % E-modulus [Pa]

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
L_BE = 0.005 %[m]
L_BD = 0.005

% Calculate I for the different intervals
x_values = linspace(0, L_ED, res);

for i = 1:res
    x = x_values(i);

    if x < (L_BE)
        d = d_E;
    elseif x < (L_BE + L_EG3 - L_BE/2 + (b_s2 / 2))
        d = d_S22;
    elseif x<(L_BE + L_EG3 - L_BE/2 + (b_s2 / 2) + L_45)
        d= d_45;
    elseif x < (L_BE + L_EG3 - L_BE/2 + (b_s2 / 2) + L_45 + b_s1 + L_G2D-(b_s1/2) - L_BD/2 )
        d=d_S21;
    else
        d = d_D;
    end

    I_shaft(i) = (pi * d^4) / 64;
end

EI = I_shaft * E;


% Integration between L_E and L_ED
for i = 2:res
    dx = x_values(i) - x_values(i-1);

    % Integrate to find rotation (omega)
    theta(i) = theta(i-1) + (M(i) / EI(i)) * dx;

    % Integrate to find deflection
    delta(i) = delta(i-1) + theta(i) * dx;
end

% Apply boundary conditions (deflection at bearings is zero), deflection is 0 at bearing E and D
index_L_ED = res;

% Calculate the correction factor K_3
K_3 = delta(index_L_ED) / L_ED;

% Correct the deflection
delta_corrected = delta - K_3 * x_values;
theta_corrected = theta - K_3;

maxDeflection = max( abs(delta_corrected) );
checkEmpiricalRequirement = maxDeflection / L_ED;

if checkEmpiricalRequirement <= 1/3000
    disp("Deflection Good")
else
    disp("Deflection not good")
end

maxTheta = max(theta_corrected);
checkThetaRequirement = tan(maxTheta) ;

if checkThetaRequirement < 0.001
    disp("Theta Good")
else
    disp("Theta Not Good")
end

% Plot the results
figure;
hold on;
delta1 = plot(x_values, delta, '--r', 'LineWidth', lwDeflection);
delta2 = plot(x_values, delta_corrected, 'r', 'LineWidth', lwDeflection);
plot(x_values(1), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(res), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
xlabel('Length [m]')
ylabel('Deflection [m]')
title('\textbf{Deflection $\delta$ of shaft 2}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend([delta1, delta2], 'No Correction', 'Corrected', ...
            'location', 'northwest')
grid on;

% Convert to degrees
theta = theta * 180/pi; % [degrees]
theta_corrected = theta_corrected * 180/pi; % [degrees]

figure;
hold on;
plot(x_values, theta, '--k', 'DisplayName', 'No correction', ...
    'LineWidth', lwDeflection)
plot(x_values, theta_corrected, 'k', 'DisplayName', 'Corrected', ...
    'LineWidth', lwDeflection)
xlabel('Length [m]')
ylabel('Angle [degrees]')
title('\textbf{Beam Slope $\theta$ of shaft 2}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend('location', 'northwest')
grid on;