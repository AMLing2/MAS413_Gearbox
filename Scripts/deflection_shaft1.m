clc; clear; close all;

% Import from Shaft 1
load('loadingDiagram_shaft1.mat')

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

% Calculate I for the different intervals
x_values = linspace(0, L_AC, res);

for i = 1:res
    x = x_values(i);
    if     x < (L_AB + (b_B/2) )
        d = d_B;
    elseif x < ( L_AG1 + (b_s1/2) )
        d = d_S1;
    elseif x < ( L_AG1 + (b_s1/2) + L_12 )
        d = d_12;
    else
        d = d_C;
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
plot(x_values(index_L_AB), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(index_L_AC), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
xlabel('Length [m]')
ylabel('Deflection [m]')
title('\textbf{Deflection $\delta$ of shaft 1}', 'interpreter', ...
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
title('\textbf{Beam Slope $\theta$ of shaft 1}', 'interpreter', ...
        'latex', 'FontSize', sizeDeflectionText)
legend('location', 'northwest')

sgtitle('Forced Deflection of Shaft 1')

pos = get(fig, 'Position'); % [xPos yPos w h]
pos(3) = pos(3)*2;
pos(1) = pos(1)*2/5;
set(fig, 'Position', pos); % set double width of default

%% Critical speed calculations

index_L_AG1 = find(x_values >= L_AG1, 1, 'first');

%abs to make sure that imaginary numbers dont get calculated in omega_c
delta_g_G1 = abs(delta_G_corrected(index_L_AG1));

% Machine Design - Equation (10.25c) page 636
omega_c = sqrt(g * (   ( (mass_g1*delta_g_G1) / (mass_g1*delta_g_G1^2))   ));
n_c = (60/(2*pi))* omega_c; % [RPM]

disp('===== Lateral Deflection =====')

% Test if n_shaft1 is outside [0.8*n_c, 1.25*n_c]
if (n_1 < 0.8 * n_c) || (n_1 > 1.25 * n_c)
    disp("Lateral vibration good");
else
    disp("Lateral vibration not good, adjusted diameter values");
    k = sqrt(  n_1*4   /    n_c     );


    d_B = d_B       *k;
    d_S1 = d_S1     *k;
    d_12 = d_12     *k;
    d_C = d_C       *k;
    

end



figure
hold on; grid on
plot(x_values, delta_G_corrected, 'r', 'LineWidth', lwDeflection);
plot(x_values(index_L_AB), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
plot(x_values(index_L_AC), 0, 'ok', 'MarkerSize', ok, 'LineWidth',1.2)
title('Lateral Deflection of Shaft 1')
xlabel('Length [m]')
ylabel('Deflection [m]')


save("deflection_shaft1.mat","d_C","d_12","d_S1","d_B")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formula for calculating new diameter based on ciritcal speed%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Needs to be done for every diameter of the shaft

%n_c_new = n_1*4;

%n_c_old = n_c

%d_new = d_old * sqrt(  n_c_new   /    n_c_old     );