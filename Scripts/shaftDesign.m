close all; clear; clc;

% TO DO:
% - Create function for diameter equation with option for 1st itteration
% - Lecture 2 slide 6-8, how to implement?
% - Modified-Godman
%   * Restructure as function in separate file
%   * Adjust formatting to match other plots
%   * Move x- and y-labels to the ends of axis
%   * Fill inn missing for safety factors


% Input parameters
n_f = 2; % Safety factor
material = 1045; % (1045 4130 4140 4340)    ! PLACEHOLDER
d_shaft = 167; % [mm] Shaft diameter        ! PLACEHLDER VALUE
r_fillet = 1; % [mm] Fillet radius          ! PLACEHLDER VALUE
D_d = 1.2; %                                ! PLACEHLDER VALUE
load_type = "Complex axial";  % ("Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial");
surface_finish = "Machined"; % ("Ground" "Machined" "Hot-rolled" "As-forged") For other types, see Machine Design page 368, Figure 6-26
reliability = 99; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)
operating_temperature = 70; % Celsius, defined by Kjell (only significant if > 450)

S_yc = -800; % [Mpa] Complressive yeild strength ! PLACEHOLDER VALUE must be incorporated with material data table

% Calculated values
R = (d/2); % [mm] Shaft radius
A = pi*R^2; % [mm^2] Shaft area
I = (pi/4)*R^4; % [mm^4] Moment of inertia
I_p = (pi/2)*R^4; % [mm^4] Polar moment of inertia

% Conversion factors
Mpa_to_ksi = 0.1450377377; % Mpa to ksi conversion factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From loadingDiagrams_shaft1.m  ! MUST BE UPDATED 

% Axial
P_x = -1; % [N]

% Shear
V_xy = 2; % [N]
V_xz = 90; % [N]

% Bending
M_z = 1*1e3; % [Nmm]
M_y = 13*1e3; % [Nmm]

% Tourqe
T = 83; % [Nmm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Axial (constant)
P_x_max = P_x; % [N]
P_x_min = P_x; % [N]

% Shear (fully reversed)
V_xy_max =  V_xy; % [N]
V_xy_min = -V_xy; % [N]
V_xz_max =  V_xz; % [N]
V_xz_min = -V_xz; % [N]

% Bending (fully reversed)
M_z_max =  M_z; % [Nmm]
M_z_min = -M_z; % [Nmm]
M_y_max =  M_y; % [Nmm]
M_y_min = -M_y; % [Nmm]

% Tourqe (constant)
T_max = T; % [Nmm]
T_min = T; % [Nmm]

%%%%% Mean & Amplitude nominal stress %%%%% (Maskinelementer, lecture 3 slide 15-18)
% Axial
sigma_x_axial_max_nom = P_x_max/A; % [Mpa]
sigma_x_axial_min_nom = P_x_min/A; % [Mpa]
sigma_x_axial_mean_nom = (sigma_x_axial_max_nom + sigma_x_axial_min_nom)/2; % [Mpa]
sigma_x_axial_amp_nom =  (sigma_x_axial_max_nom - sigma_x_axial_min_nom)/2; % [Mpa]

% Shear
tau_xy_shear_max_nom = (4/3)*(V_xy_max/A); % [Mpa]
tau_xy_shear_min_nom = (4/3)*(V_xy_min/A); % [Mpa]
tau_xy_shear_mean_nom = (tau_xy_shear_max_nom + tau_xy_shear_min_nom)/2; % [Mpa]
tau_xy_shear_amp_nom =  (tau_xy_shear_max_nom - tau_xy_shear_min_nom)/2; % [Mpa]

tau_xz_shear_max_nom = (4/3)*(V_xz_max/A); % [Mpa]
tau_xz_shear_min_nom = (4/3)*(V_xz_min/A); % [Mpa]
tau_xz_shear_mean_nom = (tau_xz_shear_max_nom + tau_xz_shear_min_nom)/2; % [Mpa]
tau_xz_shear_amp_nom =  (tau_xz_shear_max_nom - tau_xz_shear_min_nom)/2; % [Mpa]

% Bending
sigma_x_bend_max_nom = (M_max*R)/I; % [Mpa]
sigma_x_bend_min_nom = (M_min*R)/I; % [Mpa]
sigma_x_bend_mean_nom = (sigma_x_bend_max_nom + sigma_x_bend_min_nom)/2; % [Mpa]
sigma_x_bend_amp_nom =  (sigma_x_bend_max_nom - sigma_x_bend_min_nom)/2; % [Mpa]

% Torsion
tau_xy_tor_max_nom = (T_max*R)/I_p; % [Mpa]
tau_xy_tor_min_nom = (T_min*R)/I_p; % [Mpa]
tau_xy_tor_mean_nom = (tau_xy_tor_max_nom + tau_xy_tor_min_nom)/2; % [Mpa]
tau_xy_tor_amp_nom =  (tau_xy_tor_max_nom - tau_xy_tor_min_nom)/2; % [Mpa]

tau_xz_tor_max_nom = (T_max*R)/I_p; % [Mpa]
tau_xz_tor_min_nom = (T_min*R)/I_p; % [Mpa]
tau_xz_tor_mean_nom = (tau_xz_tor_max_nom + tau_xz_tor_min_nom)/2; % [Mpa]
tau_xz_tor_amp_nom =  (tau_xz_tor_max_nom - tau_xz_tor_min_nom)/2; % [Mpa]


%%%%% Material data %%%%% (Machine Design, Table A8 & A9 page 1039-1040)
material_key = [1045, 4130, 4140, 4340];
material_data = struct('S_y', num2cell([593, 655, 1462, 1365]),...
                       'S_ut', num2cell([779, 862, 1627, 1627]));
material_table = dictionary(material_key, material_data);
S_y = material_table(material).S_y;
S_ut = material_table(material).S_ut;


%%%%% Neubler's Constant for Steels %%%%% (Machine Design, Table 6-6 page 382)
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240];% [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi_temp = S_ut * Mpa_to_ksi; % [ksi] converts S_ut from Mpa to ksi
S_ut_ksi = find_closest_value(S_ut_ksi_temp, S_ut_ksi_table);

a_sqrt_in = Neublers_table(S_ut_ksi); % [in^1/2]
a_sqrt_mm = a_sqrt_in * sqrt(25.4); % [mm^1/2] Usikker på om dette er rett. (Maskinelementer, lecture 3 slide 13) 

% Notch sensitivity factor (Maskinelementer, lecture 3 slide 13 & Machine Design, equation 6.13 page 381 )
q = 1/(1 + (a_sqrt_mm/sqrt(r_fillet))); 


%%%%% Stress concentration factors %%%%% (Maskinelementer, lecture 3 slide 12)
% Geometrical (theoretical) stress concentration factors:
% See appendix C for tables and values (side 1048-1049)
D_d_bend_key = [6.00, 3.00, 2.00, 1.50, 1.20, 1.10, 1.07, 1.05, 1.03, 1.02, 1.01];
D_d_bend_values = struct('A', num2cell([0.87868, 0.89334, 0.90879, 0.93836, 0.97098, 0.95120, 0.97527, 0.98137, 0.98061, 0.96048, 0.91938]),...
                         'b', num2cell([-0.33243, -0.30860, -0.28598, -0.25759, -0.21796, -0.23757, -0.20958, -0.19653, -0.18381, -0.17711, -0.17032]));
D_d_bend_table = dictionary(D_d_bend_key, D_d_bend_values);

A_bend = D_d_bend_table(D_d).A;
b_bend = D_d_bend_table(D_d).b;

D_d_tor_key = [2.00, 1.33, 1.20, 1.09];
D_d_tor_values = struct('A', num2cell([0.86331, 0.84897, 0.83425, 0.90337]),...
                        'b', num2cell([-0.23865, -0.23161, -0.21649, -0.12692]));
D_d_tor_table = dictionary(D_d_tor_key, D_d_tor_values);

A_tor = D_d_tor_table(D_d).A;
b_tor = D_d_tor_table(D_d).b;

D_d_axial_key = [2.00, 1.50, 1.30, 1.20, 1.15, 1.10, 1.07, 1.05, 1.02, 1.01];
D_d_axial_values = struct('A', num2cell([1.01470, 0.99957, 0.99682, 0.96272, 0.98084, 0.98450, 0.98498, 1.00480, 1.01220, 0.98413]), ...
                          'b', num2cell([-0.30035, -0.28221, -0.25751, -0.25527, -0.22485, -0.20818, -0.19548, -0.17076, -0.12474, -0.10474]));
D_d_axial_table = dictionary(D_d_axial_key, D_d_axial_values);

A_axial = D_d_axial_table(D_d).A;
b_axial = D_d_axial_table(D_d).b;

K_t_bend = A_bend*(r_fillet/d_shaft)^b_bend;  % for normal stress, Appendix C
K_t_tor = A_tor*(r_fillet/d_shaft)^b_tor; % for shear stress
K_t_axial = A_axial*(r_fillet/d_shaft)^b_axial;

% Fatigue (dynamic) stress concentration foactors:
K_f_bend = 1 + q * (K_t_bend - 1);    % for normal stress, Appendix C
K_f_tor = 1 + q * (K_t_tor - 1);    % for shear stress
K_f_axial = 1 + q * (K_t_axial - 1); %  for ??


%%%%% Mean & Amplitude stress with stress concentration (Ductile materials) %%%%% (Maskinelementer, lecture 3 slide 19 & 21)
% Axial
sigma_x_axial_mean = sigma_x_axial_mean_nom * K_f_axial; % [Mpa]
sigma_x_axial_amp =  sigma_x_axial_amp_nom  * K_f_axial; % [Mpa]

% Shear
tau_xy_shear_mean = tau_xy_shear_mean_nom * K_f_tor; % [Mpa]
tau_xy_shear_amp =  tau_xy_shear_amp_nom  * K_f_tor; % [Mpa]
tau_xz_shear_mean = tau_xz_shear_mean_nom * K_f_tor; % [Mpa]
tau_xz_shear_amp =  tau_xz_shear_amp_nom  * K_f_tor; % [Mpa]

% Bending
sigma_x_bend_mean = sigma_x_bend_mean_nom * K_f_bend; % [Mpa]
sigma_x_bend_amp =  sigma_x_bend_amp_nom  * K_f_bend; % [Mpa]

% Torsion
tau_xy_tor_mean = tau_xy_tor_mean_nom * K_f_tor; % [Mpa]
tau_xy_tor_amp =  tau_xy_tor_amp_nom  * K_f_tor; % [Mpa]
tau_xz_tor_mean = tau_xz_tor_mean_nom * K_f_tor; % [Mpa]
tau_xz_tor_amp =  tau_xz_tor_amp_nom  * K_f_tor; % [Mpa]

% Resultant mean and amplitude (Maskinelementer, lecture 3 slide 23-24)
simga_x_mean = sigma_x_axial_mean + sigma_x_bend_mean; % [Mpa]
simga_x_amp =  sigma_x_axial_amp + sigma_x_bend_amp;   % [Mpa]

tau_xy_mean = tau_xy_shear_mean + tau_xy_tor_mean; % [Mpa]
tau_xy_amp =  tau_xy_shear_amp + tau_xy_tor_amp;   % [Mpa]
tau_xz_mean = tau_xz_shear_mean + tau_xz_tor_mean; % [Mpa]
tau_xz_amp =  tau_xz_shear_amp + tau_xz_tor_amp;   % [Mpa]

% Von Mises for
% sigma_von_mises = sqrt((sigma_x-sigma_y)^2 + (sigma_y-sigma_z)^2 + (sigma_z-sigma_x)^2 + 6*(tau_xy^2+tau_yz^2tau_zx^2)/2);

%%%%% Correction factors %%%%%
% Load factor % (Maskinelementer, lecture 3 slide 43 & Machine Design, page 366)
C_load_table_key = ["Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial"];
C_load_table_value = [1 0.75 1 0.75 1];
C_load_table = dictionary(C_load_table_key, C_load_table_value);

C_load = C_load_table(load_type);

% Size factor (Maskinelementer, lecture 3 slide 44 & Machine Design, page 367)
if d_shaft <= 8
    C_size = 1;
elseif d_shaft <= 250
    C_size = 1.189*d_shaft^(-0.097);
else
    C_size = 0.6;
end

% Surface factor (Maskinelementer, lecture 3 slide 45 & Machine Design, page 369)
C_surf_table_key = ["Ground" "Machined" "Hot-rolled" "As-forged"];
C_surf_table_value = struct('A', num2cell([1.58 4.51 57.7 272]), ...
                            'b', num2cell([-0.085 -0.265 -0.718 -0.995]));
C_surf_table = dictionary(C_surf_table_key, C_surf_table_value);
C_surf_A = C_surf_table(surface_finish).A;
C_surf_b = C_surf_table(surface_finish).b;

C_surf = C_surf_A*S_ut^C_surf_b;
if C_surf > 1
    C_surf = 1;
end

% Temperature factor (Maskinelementer, lecture 3 slide 47 & Machine Design, page 371)
if operating_temperature <= 450
    C_temp = 1;
elseif operating_temperature <= 550
    C_temp = 1-0.0058*(operating_temperature-450);
else
    fprintf("Error: Operating temperature too high\n")
end

% Reliability factor (Maskinelementer, lecture 3 slide 48 & Machine Design, page 371)
C_reliab_table_key = [50 90 95 99 99.9 99.99 99.999 99.9999];
C_reliab_table_value = [1 0.897 0.868 0.814 0.753 0.702 0.659 0.620];
C_reliab_table = dictionary(C_reliab_table_key, C_reliab_table_value);
C_reliab = C_reliab_table(reliability);

% Other factors
% See lecture 3 slide 49

% For steels with "knee" (Maskinelementer, lecture 3 slide 39)
if S_ut < 1400
    S_e_prime = 0.5*S_ut;
else
    S_e_prime = 700;
end

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % (Maskinelementer, lecture 4 slide 5)


% Diameter eq......
d = ((16*n_f/pi)*(sqrt(4*(K_f_bend*M_amp)^2+3*(K_f_tor*T_amp)^2)/S_e)+(sqrt(4*(K_f_bend*M_mean)^2+3*(K_f_tor*T_mean)^2/S_ut)))^(1/3)

% Quick check: failure againt yels at the first cycle (Maskinelementer, lecture 5 slide 7) 
sigma_prime_amp = sqrt(((32*K_f_bend*M_amp)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_amp)/(pi*d_shaft^3)));
sigma_prime_mean = sqrt(((32*K_f_bend*M_mean)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_mean)/(pi*d_shaft^3))); 
sigma_max = sigma_prime_mean + sigma_prime_amp;
n_y = S_y / sigma_max

%% First itteration

% Estimating stress geometric concentration factors for preliminary stage (Maskinelementerlecture 5 slide 10)
% Shoulder fillet sharp (r/d = 0.02, D/d = 1.5)
% K_t_bend = 2.7;
% K_t_tor = 2.2;
% K_t_axial = 3.0;

% Shoulder fillet well-rounded (r/d = 0.1, D/d = 1.5)
K_t_bend = 1.7;
K_t_tor = 1.5;
K_t_bend = 1.9;

% End-mill keyset (r/d = 0.02)
% K_t_bend = 2.14;
% K_t_tor = 3;
% K_t_bend = -; 

% Sled runner keyset
% K_t_bend = 1.7;
% K_t_tor = -;
% K_t_bend = -;

% Retaining ring groove
% K_t_bend = 5;
% K_t_tor = 3;
% K_t_bend = 5;

% Conservative estimate for preliminary stage (q is unknown)
K_f_bend = K_t_bend;
K_f_tor = K_t_tor;
K_f_axial = K_t_bend;

% Correction factors for preliminary stage % (Maskinelementer, lecture 5 slide 12)
C_load = 1; % Bending
C_size = 1;
% C_surf = C_surf;
% C_temp = C_temp; 
% C_reliab = C_reliab;

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % (Maskinelementer, lecture 4 slide 5)

% Shaft diameter
d1 = ((16*n_f/pi)*(sqrt(4*(K_f_bend*M_amp)^2+3*(K_f_tor*T_amp)^2)/S_e)+(sqrt(4*(K_f_bend*M_mean)^2+3*(K_f_tor*T_mean)^2/S_ut)))^(1/3)


%% Modified-Goodman Graph (Work in progress)
% (Maskinelementer, lecture 4 slide 12-22)
close all;

% Define mean stress range
sigma_mean_goodman_L = linspace(S_yc, 0, 100); % Left side (compressive)
sigma_mean_goodman_R = linspace(0, S_ut, 100); % Right side (tensile)

% Modified-Goodman equation
sigma_amp_goodman_L = S_e * ones(size(sigma_mean_goodman_L));
sigma_amp_goodman_R = S_e * (1 - sigma_mean_goodman_R / S_ut);

% Static yeilding line (first cycle)
S_y_goodman_L = S_y * (1 - sigma_mean_goodman_L / S_yc);
S_y_goodman_R = S_y * (1 - sigma_mean_goodman_R / S_y);

% Check for static failure
n_y = S_y/simga_max;

% Check for fatigue failure
sigma_rev = sigma_amp_vm/(1-(sigma_mean_vm/S_ut));
n_f = S_e/simga_rev;

% Intersecting points for indexing
[~, intersect_L] = min(abs(sigma_amp_goodman_L - S_y_goodman_L));
[~, intersect_R] = min(abs(sigma_amp_goodman_R - S_y_goodman_R));

figure; hold on;
% Colour area between graphs and x-axis
fill([sigma_mean_goodman_L, sigma_mean_goodman_R], [S_y_goodman_L(1:intersect_L-1), sigma_amp_goodman_L(intersect_L:end),...
 sigma_amp_goodman_R(1:intersect_R-1), S_y_goodman_R(intersect_R:end)], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;

% Equations
plot(sigma_mean_goodman_L, sigma_amp_goodman_L, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2); hold on
plot(sigma_mean_goodman_R, sigma_amp_goodman_R, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
plot(sigma_mean_goodman_L, S_y_goodman_L, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
plot(sigma_mean_goodman_R, S_y_goodman_R, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)

% Points for strengths
plot(0, S_e, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_e point
plot(S_ut, 0, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_ut point
plot(0, S_y, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
plot(S_y, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
plot(S_yc, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_yc point
axis([S_yc-100 S_ut+100 0 S_y+100]);

% Axis labels ! NEEDS ADJUSTMENT
xticks([S_yc, S_y, S_ut]); % Set x-axis tick positions
yticks([S_e, S_y]); % Set y-axis tick positions
xticklabels({'S_{yc}', 'S_{y}', 'S_{ut}'}); % Custom x-axis labels
yticklabels({'S_e', 'S_y'}); % Custom y-axis labels
ax = gca;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ax.YAxisLocation = 'origin'; % Move y-axis to x = 0

xlabel('\sigma_m [MPa]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma_a [MPa]', 'FontSize', 14, 'FontWeight', 'bold');


%% Notes

% sigma_x_nom = sigma_x_axial+sigma_x_bending;
% 
% sigma_y_axial_nom = abs(P)/A;
% sigma_y_bend_nom = abs(M*R)/I;
% sigma_y_nom = sigma_x_axial+sigma_x_bending;
% 
% sigma_z_axial_nom = abs(P)/A;
% sigma_z_bend_nom = abs(M*R)/I;
% sigma_z_nom = sigma_x_axial+sigma_x_bending;

% Stresses with stress concentration factor
% sigma_x_axial = sigma_x_axial_nom * K_f_axial;
% sigma_x_bend = sigma_x_bend_nom * K_f_bend;
% sigma_x = sigma_x_axial + sigma_x_bend;
% 
% sigma_y_axial = sigma_y_axial_nom * K_f_axial;
% sigma_y_bend = sigma_y_bend_nom * K_f_bend;
% sigma_y = sigma_y_axial + sigma_y_bend;
% 
% sigma_z_axial = sigma_z_axial_nom * K_f_axial;
% sigma_z_bend = sigma_z_bend_nom * K_f_bend;
% sigma_z = sigma_z_axial + sigma_z_bend;
% 
% tau_shear = tau_shear_nom * K_f_tor;
% tau_tor = tau_tor_nom * K_f_tor;

% Von Mises for state of stress (Maskinelementer, lecture 2 slide 32)
% sigma_von_mises = sqrt((sigma_x-sigma_y)^2 + (sigma_y-sigma_z)^2 + (sigma_z-sigma_x)^2 + 6*(tau_xy^2+tau_yz^2tau_zx^2)/2);


% (Maskinelementer, lecture 5 slide 3)
% sigma_amp = K_f-bend*(32*M_amp)/(pi*d_shaft^3);
% sigma_mean = K_f-bend*(32*M_mean)/(pi*d_shaft^3);

% tau_amp = K_f-tor*(16*T_amp)/(pi*d_shaft^3);
% tau_mean = K_f-tor*(16*T_mean)/(pi*d_shaft^3);

% sigma_max = ?
% sigma_min = ?
% sigma_mean = (sigma_max + sigma_min)/2; % Mean (midrange) stess
% if sigma_mean ~= 0
%     fprintf("sigma_mean = ", sigma_mean, " Non nully reversed loading!")
% end
% sigma_amp = (sigma_max - sigma_min)/2; % Alternating stress amplitude


% (Maskinelementer, lecture 5 slide 6)
% sigma_rev = sigma_amp_prime/(1-(simga_mean_prime/S_ut)); % 
% n_f = S_e / sigma_rev; %
%%



% Function to find the closest value rounded down in the list (Most conservative approach)
function closest_value = find_closest_value(value, list)
    list = list(list <= value); % Filter out values greater than the given value
    if isempty(list)
        error('No values in the list are less than or equal to the given value.');
    end
    [~, index] = min(abs(list - value)); % Find the closest value among the remaining ones
    closest_value = list(index);
end