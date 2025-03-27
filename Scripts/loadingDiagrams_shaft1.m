clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAS413 Project: Loading Diagrams - Shaft 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants

g = 9.81; %[m/s^2]
m_G1 = 0.52675; %[kg]
E = 210e9; % E-modulus [Pa]

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
L_AB  = 0.05; % [m]

% Bearing widths
b_B = 30e-3; % [m] catalogue circa 16 - 47 [mm] <-- WIP
b_C = b_B; % [m]

% Import from Gear Sizing
load('gear_sizes.mat', 'd_g1', 'd_g2', 'd_g3', 'd_g4', ...
    'b_s1', 'b_s2', 'i_tot')
    
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
F_G1 = m_G1*g;
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
xz_x = [xz_x, x];
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, zeros(size(x))]; % [N]
xz_M = [xz_M, zeros(size(x))]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]
xz_Mg =[xz_Mg, zeros(size(x))];


% L_AB < x < L_AG1
x = linspace(L_AB, L_AG1, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, zeros(size(x))]; % [N]
xz_V = [xz_V, ones(size(x)) * (-F_Bz)]; % [N]
xz_M = [xz_M, ( - F_Bz*(x - L_AB) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * T_M]; % [Nm]
xz_Mg =[xz_Mg, +( F_Bg1 * (x-L_AB)) ]; %[Nm]


% L_AG1 < x < L_AC
x = linspace(L_AG1, L_AC, resolution);
xz_x = [xz_x, x];
xz_P = [xz_P, ones(size(x)) * (-F_a1)]; % [N]
xz_V = [xz_V, ones(size(x)) * (F_r1 - F_Bz)]; % [N]
xz_M = [xz_M, ( F_r1*(x - L_AG1) - F_Bz*(x - L_AB) - F_a1*(r_G1) )]; % [Nm]
xz_T = [xz_T, ones(size(x)) * (T_M - F_t1*r_G1)]; % [Nm]
xz_Mg =[xz_Mg, +( F_Bg1 * (x - L_AB) ) - ( F_G1 * (x - L_AG1) ) ]; %[Nm]


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

% Cross Section A
[~, cs_A_idx] = closest(xy_x, 0);
cs_A_P = xy_P(cs_A_idx);
cs_A_T = xy_T(cs_A_idx);
cs_A_M = M(cs_A_idx);
cs_A = [cs_A_P cs_A_T*1e3 cs_A_M*1e3]; % [(N) (Nmm) (Nmm)]

% Cross Section B
[~, cs_0_idx] = closest(xy_x, L_AB);
cs_0_P = xy_P(cs_0_idx);
cs_0_T = xy_T(cs_0_idx);
cs_0_M = M(cs_0_idx);
cs_0 = [cs_0_P cs_0_T*1e3 cs_0_M*1e3]; % [(N) (Nmm) (Nmm)]

% Cross Section 1
[~, cs_1_idx] = closest(xy_x, L_A1);
cs_1_P = xy_P(cs_1_idx);
cs_1_T = xy_T(cs_1_idx);
cs_1_M = M(cs_1_idx);
cs_1 = [cs_1_P cs_1_T*1e3 cs_1_M*1e3]; % [(N) (Nmm) (Nmm)]

% Cross Section 2
[~, cs_2_idx] = closest(xy_x, L_A2);
cs_2_P = xy_P(cs_2_idx);
cs_2_T = xy_T(cs_2_idx);
cs_2_M = M(cs_2_idx);
cs_2 = [cs_2_P cs_2_T*1e3 cs_2_M*1e3]; % [(N) (Nmm) (Nmm)]

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


%% Shaft Deflection - Forced and free
% Visuals
lwDeflection = 2;
sizeDeflectionText = 16;
ok = 10;

res = 300;

% Initialization
theta = zeros(1, res);
theta_G = zeros(1, res);
delta = zeros(1, res);
delta_G = zeros(1, res);
I_shaft = zeros(1, res);

% Diameters of shaft
d_c   = 0.010; % [m]
d_12  = 0.015; % [m]
d_S11 = 0.011; % [m]
d_S12 = 0.013; % [m]

% Calculate I for the different intervals
x_values = linspace(0, L_AC, res);

for i = 1:res
    x = x_values(i);

    if     x < (L_AB + (b_B/2) )
        d = d_S11;
    elseif x < ( L_AG1 + (b_s1/2) )
        d = d_S12;
    elseif x < ( L_AG1 + (b_s1/2) + L_12 )
        d = d_12;
    else
        d = d_c;
    end

    I_shaft(i) = (pi * d^4) / 64;
end

EI = I_shaft * E;


% Integration between L_AB and L_AC
for i = 2:res
    dx = x_values(i) - x_values(i-1);

    % Integrate to find rotation (omega)
    theta(i) = theta(i-1) + (M(i) / EI(i)) * dx;
    theta_G(i) = theta_G(i-1) + (xz_Mg(i) / EI(i)) * dx;

    % Integrate to find deflection
    delta(i) = delta(i-1) + theta(i) * dx;
    delta_G(i) = delta_G(i-1) + theta_G(i) * dx;
    
end

% Apply boundary conditions (deflection at bearings is zero), deflection is 0 at L_AB and L_AC
index_L_AB = find(x_values >= L_AB, 1, 'first');
index_L_AC = res;

% Correction Factor K_4: no deflection @ first bearing
K_4 = 0;
K_4_G = 0;

% Correction Factor K_3: non deflection @ second bearing
K_3 = delta(index_L_AC) / L_BC;
K_3_G = delta_G(index_L_AC) / L_BC; 

% Correct the Deflection and Beam Slope
delta_corrected = delta - K_3 * (x_values - x_values(index_L_AB)) - K_4;
delta_G_corrected = delta_G - K_3_G * (x_values - x_values(index_L_AB)) - K_4;
theta_corrected = theta - K_3;
theta_G_corrected = theta_G -K_3_G;

maxDeflection = max( abs(delta_corrected) );
checkEmpiricalRequirement = maxDeflection / L_AC;

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
plot(x_values(index_L_AB), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(index_L_AC), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
xlabel('Length [m]')
ylabel('Deflection [m]')
title('\textbf{Deflection $\delta$ of shaft 1}', 'interpreter', ...
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
title('\textbf{Beam Slope $\theta$ of shaft 1}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend('location', 'northwest')
grid on;

%% Critical speed calculations
%masses of the gears

index_L_AG1 = find(x_values >= L_AG1, 1, 'first');

delta_g_G1 = abs(delta_G_corrected(index_L_AG1)); %abs to make sure that imaginary numbers dont get calculated in omega_c

omega_c = sqrt(g * (   ( (m_G1*delta_g_G1) / (m_G1*delta_g_G1^2))   )); %Machine design equation 10.25c
n_c = (60/(2*pi))* omega_c; %[rpm]

% Test if n_shaft1 is outside [0.8*n_c, 1.25*n_c]
if (n_1 < 0.8 * n_c) || (n_1 > 1.25 * n_c)
    disp("Lateral vibration good");
else
    disp("Lateral vibration not good");
end

close all
figure
plot(x_values,delta_G_corrected)