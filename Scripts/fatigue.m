% Conversion factors
MPa_to_ksi = 0.1450377377; % Mpa to ksi conversion factor

% Calculated values
d_shaft = cs(1); % [mm] Shaft diameter
r_fillet = cs(2);
D_d = cs(3);
R = (d_shaft/2); % [mm] Shaft radius
A = pi*R^2; % [mm^2] Shaft area
I = (pi/4)*R^4; % [mm^4] Moment of inertia
I_p = (pi/2)*R^4; % [mm^4] Polar moment of inertia

% Forces
P_x = cs(4); % [N] % Axial
M_z = cs(5); % [Nmm] % Bending
M_y = cs(6); % [Nmm] % Bending
T = cs(7); % [Nmm] % Tourqe

V_xy = 0; % [N] % Shear
V_xz = 0; % [N] % Shear

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

M_max =  sqrt(M_z_max^2 + M_y_max^2); % [Nmm]
M_min = -sqrt(M_z_min^2 + M_y_min^2); % [Nmm]
M_mean = (M_max + M_min)/2; % [Nmm]
M_amp =  (M_max - M_min)/2; % [Nmm]

% Tourqe (constant)
T_max = T; % [Nmm]
T_min = T; % [Nmm]
T_mean = (T_max + T_min)/2; % [Nmm]
T_amp =  (T_max - T_min)/2; % [Nmm]


%%%%% Mean & Amplitude nominal stress %%%%% (Maskinelementer, lecture 3 slide 15-18)
% Axial
sigma_x_axial_max_nom = abs(P_x_max)/A; % [Mpa]
sigma_x_axial_min_nom = abs(P_x_max)/A; % [Mpa]
sigma_x_axial_mean_nom = (sigma_x_axial_max_nom + sigma_x_axial_min_nom)/2; % [Mpa]
sigma_x_axial_amp_nom =  (sigma_x_axial_max_nom - sigma_x_axial_min_nom)/2; % [Mpa]

% Shear
% tau_xy_shear_max_nom = (4/3)*(abs(V_xy_max)/A); % [Mpa]
% tau_xy_shear_min_nom = (4/3)*(V_xy_min/A); % [Mpa]
% tau_xy_shear_mean_nom = (tau_xy_shear_max_nom + tau_xy_shear_min_nom)/2; % [Mpa]
% tau_xy_shear_amp_nom =  (tau_xy_shear_max_nom - tau_xy_shear_min_nom)/2; % [Mpa]

% tau_xz_shear_max_nom = (4/3)*(abs(V_xz_max)/A); % [Mpa]
% tau_xz_shear_min_nom = (4/3)*(V_xz_min/A); % [Mpa]
% tau_xz_shear_mean_nom = (tau_xz_shear_max_nom + tau_xz_shear_min_nom)/2; % [Mpa]
% tau_xz_shear_amp_nom =  (tau_xz_shear_max_nom - tau_xz_shear_min_nom)/2; % [Mpa]

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
material_data = struct('S_y', num2cell([0, 0, 655, 0]),...
                       'S_ut', num2cell([0, 0, 1020, 0]));
material_table = dictionary(material_key, material_data);
S_y = material_table(material).S_y;
S_ut = material_table(material).S_ut;
S_yc = -S_y; % [Mpa] Compressive yeild strength ! PLACEHOLDER VALUE must be incorporated with material data table


%%%%% Neubler's Constant for Steels %%%%% (Machine Design, Table 6-6 page 382)
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240];% [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi_temp = S_ut * MPa_to_ksi; % [ksi] converts S_ut from Mpa to ksi
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

% For end-milled keyseat (Machine Design, figure 10-16, page 615)
if keyseat == "y"
    K_t_tor = K_t_keyseat;
end

% Fatigue (dynamic) stress concentration foactors:
K_f_bend = 1 + q * (K_t_bend - 1);    % for normal stress, Appendix C
K_f_tor = 1 + q * (K_t_tor - 1);    % for shear stress
K_f_axial = 1 + q * (K_t_axial - 1); %  for axial sress


%%%%% Mean & Amplitude stress with stress concentration (Ductile materials) %%%%% (Maskinelementer, lecture 3 slide 19 & 21)
% Axial
sigma_x_axial_mean = sigma_x_axial_mean_nom * K_f_axial; % [Mpa]
sigma_x_axial_amp =  sigma_x_axial_amp_nom  * K_f_axial; % [Mpa]

% Shear
% tau_xy_shear_max = tau_xy_shear_max_nom * K_f_tor; % [Mpa] Only xy !!
% tau_xy_shear_mean = tau_xy_shear_mean_nom * K_f_tor; % [Mpa]
% tau_xy_shear_amp =  tau_xy_shear_amp_nom  * K_f_tor; % [Mpa]
% tau_xz_shear_mean = tau_xz_shear_mean_nom * K_f_tor; % [Mpa]
% tau_xz_shear_amp =  tau_xz_shear_amp_nom  * K_f_tor; % [Mpa]

% Bending
% sigma_x_bend = sigma_x_bend_max_nom * K_f_bend; % [Mpa]
sigma_x_bend_mean = sigma_x_bend_mean_nom * K_f_bend; % [Mpa]
sigma_x_bend_amp =  sigma_x_bend_amp_nom  * K_f_bend; % [Mpa]

% Torsion
% tau_xy_tor_max = tau_xy_tor_max_nom * K_f_tor; % [Mpa]
tau_xy_tor_mean = tau_xy_tor_mean_nom * K_f_tor; % [Mpa]
tau_xy_tor_amp =  tau_xy_tor_amp_nom  * K_f_tor; % [Mpa]
tau_xz_tor_mean = tau_xz_tor_mean_nom * K_f_tor; % [Mpa]
tau_xz_tor_amp =  tau_xz_tor_amp_nom  * K_f_tor; % [Mpa]

% Resultant mean and amplitude (Maskinelementer, lecture 3 slide 23-24)
sigma_x_mean = sigma_x_axial_mean + sigma_x_bend_mean; % [Mpa]
sigma_x_amp =  sigma_x_axial_amp + sigma_x_bend_amp;   % [Mpa]

% Shear neglected
tau_xy_mean = tau_xy_tor_mean; % [Mpa]
tau_xy_amp =  tau_xy_tor_amp;   % [Mpa]
tau_xz_mean = tau_xz_tor_mean; % [Mpa]
tau_xz_amp =  tau_xz_tor_amp;   % [Mpa]

% Included shear (done incorrectly)
% tau_xy_mean = tau_xy_shear_mean + tau_xy_tor_mean; % [Mpa]
% tau_xy_amp =  tau_xy_shear_amp + tau_xy_tor_amp;   % [Mpa]
% tau_xz_mean = tau_xz_shear_mean + tau_xz_tor_mean; % [Mpa]
% tau_xz_amp =  tau_xz_shear_amp + tau_xz_tor_amp;   % [Mpa]

% Von Mises
sigma_vm_mean = sqrt(sigma_x_mean^2 + 3*(tau_xy_mean^2+tau_xz_mean^2)); % [Mpa]
sigma_vm_amp =  sqrt(sigma_x_amp^2 + 3*(tau_xy_mean^2+tau_xz_amp^2));   % [Mpa]
sigma_vm_max = sigma_vm_mean + sigma_vm_amp; % [Mpa]


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


%%%%% Diameter equation %%%%%
if first_itteration == "y"
    
    % Estimating stress geometric concentration factors for preliminary stage (Maskinelementerlecture 5 slide 10)
    % Shoulder fillet sharp (r/d = 0.02, D/d = 1.5)
    K_t_bend = 2.7;
    K_t_tor = 2.2;
    K_t_axial = 3.0;
    
    % Shoulder fillet well-rounded (r/d = 0.1, D/d = 1.5)
    % K_t_bend = 1.7;
    % K_t_tor = 1.5;
    % K_t_bend = 1.9;
    
    % End-mill keyset (r/d = 0.02)
    % K_t_bend = 2.14;
    % K_t_tor = 3;
    % K_t_axial = -; 
    
    % Sled runner keyset
    % K_t_bend = 1.7;
    % K_t_tor = -;
    % K_t_axial = -;
    
    % Retaining ring groove
    % K_t_bend = 5;
    % K_t_tor = 3;
    % K_t_axial = 5;
    
    % Conservative estimate for preliminary stage (q is unknown)
    K_f_bend = K_t_bend;
    K_f_tor = K_t_tor;
    K_f_axial = K_t_bend; % ??
    
    % Correction factors for preliminary stage % (Maskinelementer, lecture 5 slide 12)
    C_load = 1; % Bending
    C_size = 1;
end

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % (Maskinelementer, lecture 4 slide 5)

d_eq = ((16*n_f)/pi)*((sqrt(4*(K_f_bend*M_amp)^2+3*(K_f_tor*T_amp)^2)/S_e)+(sqrt(4*(K_f_bend*M_mean)^2+3*(K_f_tor*T_mean)^2)/S_ut))^(1/3);
fprintf('d_eq = %.2f [mm],  Recomended shaft diameter\n', d_eq)

% Quick check: failure againt yels at the first cycle (Maskinelementer, lecture 5 slide 7) 
% sigma_prime_amp = sqrt(((32*K_f_bend*M_amp)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_amp)/(pi*d_shaft^3)));
% sigma_prime_mean = sqrt(((32*K_f_bend*M_mean)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_mean)/(pi*d_shaft^3))); 
% sigma_max = sigma_prime_mean + sigma_prime_amp;
% n_y = S_y / sigma_max


%%%%% Check if shear can be neglected %%%%% INCOMPLETE
% fprintf('sigma_bend = %.2f [MPa]\n', sigma_x_bend_max_nom + sigma_x_axial_max_nom)
% fprintf('tau_tor = %.2f [MPa]\n', tau_xy_tor_max_nom)
% fprintf('tau_shear = %.2f [MPa]\n', tau_xy_shear_max_nom + tau_xz_shear_max_nom)
% shear_bending = (tau_xy_shear_max_nom + tau_xz_shear_max_nom)/(sigma_x_bend_max_nom + sigma_x_axial_max_nom);
% shear_tor = (tau_xy_shear_max_nom + tau_xz_shear_max_nom)/tau_xy_tor_max_nom;
% 
% fprintf('tau_shear / sigma_bending = %.2f --> ', shear_bending)
% if shear_bending <= 0.1
%     fprintf('shear force can be neglected\n')
% else
%     fprintf('shear force can NOT be neglected\n')
% end
% 
% fprintf('tau_shear / tau_tor = %.2f --> ', shear_tor)
% if shear_tor <= 0.1
%     fprintf('shear force can be neglected\n')
% else
%     fprintf('shear force can NOT be neglected\n')
% end



% Function to find the closest value rounded down in the list (Most conservative approach)
function closest_value = find_closest_value(value, list)
    list = list(list <= value); % Filter out values greater than the given value
    if isempty(list)
        error('No values in the list are less than or equal to the given value.');
    end
    [~, index] = min(abs(list - value)); % Find the closest value among the remaining ones
    closest_value = list(index);
end