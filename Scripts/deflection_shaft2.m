clc; clear; close all;

% Import from Shaft 2
load('loadingDiagram_shaft2.mat')

%% Shaft deflection calculations - free and forced

% Visuals
lwDeflection = 2;
sizeDeflectionText = 16;
ok = 10;

res = 300;

%Initialize arrays:
I_shaft = [];

theta = zeros(1, res);
theta_G = zeros(1, res);

theta_corrected = [];
theta_corrected_G = [];

delta = zeros(1, res);
delta_G = zeros(1, res);

delta_corrected = [];
delta_corrected_G = [];

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
    theta_G(i) = theta_G(i-1) + (xz_Mg(i) / EI(i)) * dx;

    % Integrate to find deflection
    delta(i) = delta(i-1) + theta(i) * dx;
    delta_G(i) = delta_G(i-1) + theta_G(i) * dx;
end

% Apply boundary conditions (deflection at bearings is zero), deflection is 0 at bearing E and D
index_L_ED = res;

% Calculate the correction factor K_3
K_3 = delta(index_L_ED) / L_ED;
K_3_G = delta_G(index_L_ED) / L_ED;

% Correct the deflection
delta_corrected = delta - K_3 * x_values;
delta_corrected_G = delta_G - K_3_G * x_values; 
theta_corrected = theta - K_3;
theta_corrected_G = theta_G - K_3_G;

maxDeflection = max( abs(delta_corrected) );
checkEmpiricalRequirement = maxDeflection / L_ED;

disp('===== Forced Deflection =====')

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

% Forced Deflection Plot
fig = figure;
subplot(1,2,1)
hold on; grid on
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

% Convert to degrees
theta = theta * 180/pi; % [degrees]
theta_corrected = theta_corrected * 180/pi; % [degrees]

% Beam Slope Plot
subplot(1,2,2)
hold on; grid on
plot(x_values, theta, '--k', 'DisplayName', 'No correction', ...
    'LineWidth', lwDeflection)
plot(x_values, theta_corrected, 'k', 'DisplayName', 'Corrected', ...
    'LineWidth', lwDeflection)
xlabel('Length [m]')
ylabel('Angle [degrees]')
title('\textbf{Beam Slope $\theta$ of shaft 2}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend('location', 'northwest')

sgtitle('Forced Deflection of Shaft 2')

pos = get(fig, 'Position'); % [xPos yPos w h]
pos(3) = pos(3)*2;
pos(1) = pos(1)*2/5;
set(fig, 'Position', pos); % set double width of default

%% Critical speed calculations

index_L_EG2 = find(x_values >= L_EG2, 1, 'first');
index_L_EG3 = find(x_values >= L_EG3, 1, 'first');

delta_g_G2 = abs(theta_corrected_G(index_L_EG2));
delta_g_G3 = abs(theta_corrected_G(index_L_EG3));

% Machine Design - Equation (10.25c) page 636
omega_c = sqrt(g * (  ( ( (m_G2*delta_g_G2) + (m_G3*delta_g_G3) ) / ...
                      ( (m_G2*delta_g_G2^2) + (m_G3*delta_g_G3^2) ) )));
n_c = (60/(2*pi))* omega_c; % [RPM]

n_shaft2 = n_1/i_s1;

disp('===== Lateral Deflection =====')

% Test if n_shaft2 is outside [0.8*n_c, 1.25*n_c]
if (n_shaft2 < 0.8 * n_c) || (n_shaft2 > 1.25 * n_c)
    disp("Lateral vibration good");
else
    disp("Lateral vibration not good, values adjusted");
    
    k =  sqrt(  n_shaft2*4   /    n_c     );

    d_E = d_E*k;
    d_S22= d_S22*k;
    d_45= d_45*k;
    d_S21 = d_S21*k;
    d_D = d_D*k;
    
end

save("deflection_shaft2.mat","d_D","d_S21","d_45","d_S22","d_E")


figure
hold on; grid on
plot(x_values, delta_corrected_G, 'r', 'LineWidth', lwDeflection);
plot(x_values(1), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(res), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
title('Lateral Deflection of Shaft 2')
xlabel('Length [m]')
ylabel('Deflection [m]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formula for calculating new diameter based on ciritcal speed%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Needs to be done for every diameter of the shaft

%n_c_new = n_shaft2*4;

%n_c_old = n_c

%d_new = d_old * sqrt(  n_c_new   /    n_c_old     );