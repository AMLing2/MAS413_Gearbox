clc; clear; close all;

% Import from Shaft 3
load('loadingDiagram_shaft3.mat')

%% Shaft Deflection - Free and forced

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
x_values = linspace(0, L_FH, res);

for i = 1:res
    x = x_values(i);

    if     x < ( (b_F/2) )
        d = d_F;
    elseif x < ( (b_F/2) + L_78 )
        d = d_78;
    elseif x < ( (b_F/2) + L_78 + b_s2 + L_45 + b_s1 + L_12)
        d = d_G;
    else
        d = d_S3;
    end

    I_shaft(i) = (pi * d^4) / 64;
end

EI = I_shaft * E;


% Numerical Integration
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
index_L_FG = find(x_values >= L_FG, 1, 'first');

% Correction Factor K_4: no deflection @ first bearing
K_4 = 0;

% Correction Factor K_3: non deflection @ second bearing
K_3 = delta(index_L_FG) / L_FG;
K_3_G = delta_G(index_L_FG) / L_FG;

% Correct the Deflection and Beam Slope
delta_corrected = delta - K_3 * x_values;
delta_corrected_G = delta_G - K_3_G * x_values; 
theta_corrected = theta - K_3;
theta_corrected_G = theta_G - K_3_G;

maxDeflection = max( abs(delta_corrected) );
checkEmpiricalRequirement = maxDeflection / L_FH;

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
plot(0, 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(index_L_FG), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
xlabel('Length [m]')
ylabel('Deflection [m]')
title('\textbf{Deflection $\delta$ of shaft 3}', 'interpreter', ...
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
title('\textbf{Beam Slope $\theta$ of shaft 3}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend('location', 'northwest')

sgtitle('Forced Deflection of Shaft 3')

pos = get(fig, 'Position'); % [xPos yPos w h]
pos(3) = pos(3)*2;
pos(1) = pos(1)*2/5;
set(fig, 'Position', pos); % set double width of default

%% Critical speed calculations

index_L_EG4 = find(x_values >= L_FG4, 1, 'first');

delta_g_G4 = abs(theta_corrected_G(index_L_EG4));

% Machine Design - Equation (10.25c) page 636
omega_c = sqrt( g * (  (m_G4 * delta_g_G4) / (m_G4*delta_g_G4^2)  ) );
n_c = (60/(2*pi))* omega_c; %[rpm]

n_shaft3 = n_1/i_tot;

disp('===== Lateral Deflection =====')

% Test if n_shaft1 is outside [0.8*n_c, 1.25*n_c]
if (n_shaft3 < 0.8 * n_c) || (n_shaft3 > 1.25 * n_c)
    disp("Lateral vibration good");
else
    disp("Lateral vibration not good, values adjusted");


    k =  sqrt(  n_shaft3*4   /    n_c     );
    
    d_F = d_F*k;
    d_78 = d_78*k;
    d_G * d_G*k;
    d_S3 = d_S3*k;

end


save("deflection_shaft3.mat","d_F","d_78","d_G","d_S3")



figure
hold on; grid on
plot(x_values, delta_corrected_G, 'r', 'LineWidth', lwDeflection);
plot(0, 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(index_L_FG), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
title('Lateral Deflection of Shaft 3')
xlabel('Length [m]')
ylabel('Deflection [m]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formula for calculating new diameter based on ciritcal speed if nessesary%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Needs to be done for every diameter of the shaft

%n_c_new = n_shaft3*4;

%n_c_old = n_c

%d_new = d_old * sqrt(  n_c_new   /    n_c_old     );