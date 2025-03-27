% Conversion factors
MPa_to_ksi = 0.1450377377; % [MPa] to [ksi] conversion factor

%%%%% Imported values from shaftDesign.m %%%%%
% Forces
P_x = cs(1); % [N] % Axial
T = cs(2); % [Nmm] % Torque
M = cs(3); % [Nmm] % Bending

% Neglect shear
V_xy = 0; % [N] % Shear
V_xz = 0; % [N] % Shear

% Geometry for stress concentration
d_shaft = cs(4);  % [mm] Shaft diameter
r_fillet = cs(5); % [mm] Notch fillet radius
D_d = cs(6);      % [-]
keyseat = cs(8);  % Is 1 for keyseat, 0 for shoulder-fillet

% Calculated values
R = (d_shaft/2);  % [mm] Shaft radius
A = pi*R^2;       % [mm^2] Shaft area
I = (pi/4)*R^4;   % [mm^4] Moment of inertia (I_x = I_y)
I_p = (pi/2)*R^4; % [mm^4] Polar moment of inertia


%%%%% Max / Min forces %%%%%
% Axial (constant)
P_x_max = P_x; % [N]
P_x_min = P_x; % [N]

% Shear (fully reversed)
V_xy_max =  V_xy; % [N]
V_xy_min = -V_xy; % [N]
V_xz_max =  V_xz; % [N]
V_xz_min = -V_xz; % [N]

% Bending (fully reversed)
M_max =  M; % [Nmm]
M_min = -M; % [Nmm]
M_mean = (M_max + M_min)/2; % [Nmm]
M_amp =  (M_max - M_min)/2; % [Nmm]

% Torque (constant)
T_max = T; % [Nmm]
T_min = T; % [Nmm]
T_mean = (T_max + T_min)/2; % [Nmm]
T_amp =  (T_max - T_min)/2; % [Nmm]


%%%%% Mean & Amplitude nominal stress %%%%% MAS236 L3 s15-18
% Axial
sigma_x_axial_max_nom = abs(P_x_max)/A; % [MPa]
sigma_x_axial_min_nom = abs(P_x_max)/A; % [MPa]
sigma_x_axial_mean_nom = (sigma_x_axial_max_nom + sigma_x_axial_min_nom)/2; % [MPa]
sigma_x_axial_amp_nom =  (sigma_x_axial_max_nom - sigma_x_axial_min_nom)/2; % [MPa]

% Shear
% tau_xy_shear_max_nom = (4/3)*(abs(V_xy_max)/A); % [MPa]
% tau_xy_shear_min_nom = (4/3)*(V_xy_min/A); % [MPa]
% tau_xy_shear_mean_nom = (tau_xy_shear_max_nom + tau_xy_shear_min_nom)/2; % [MPa]
% tau_xy_shear_amp_nom =  (tau_xy_shear_max_nom - tau_xy_shear_min_nom)/2; % [MPa]

% tau_xz_shear_max_nom = (4/3)*(abs(V_xz_max)/A); % [MPa]
% tau_xz_shear_min_nom = (4/3)*(V_xz_min/A); % [MPa]
% tau_xz_shear_mean_nom = (tau_xz_shear_max_nom + tau_xz_shear_min_nom)/2; % [MPa]
% tau_xz_shear_amp_nom =  (tau_xz_shear_max_nom - tau_xz_shear_min_nom)/2; % [MPa]

% Bending
sigma_x_bend_max_nom = (M_max*R)/I; % [MPa]
sigma_x_bend_min_nom = (M_min*R)/I; % [MPa]
sigma_x_bend_mean_nom = (sigma_x_bend_max_nom + sigma_x_bend_min_nom)/2; % [MPa]
sigma_x_bend_amp_nom =  (sigma_x_bend_max_nom - sigma_x_bend_min_nom)/2; % [MPa]

% Torsion
tau_xy_tor_max_nom = (T_max*R)/I_p; % [MPa]
tau_xy_tor_min_nom = (T_min*R)/I_p; % [MPa]
tau_xy_tor_mean_nom = (tau_xy_tor_max_nom + tau_xy_tor_min_nom)/2; % [MPa]
tau_xy_tor_amp_nom =  (tau_xy_tor_max_nom - tau_xy_tor_min_nom)/2; % [MPa]

tau_xz_tor_max_nom = (T_max*R)/I_p; % [MPa]
tau_xz_tor_min_nom = (T_min*R)/I_p; % [MPa]
tau_xz_tor_mean_nom = (tau_xz_tor_max_nom + tau_xz_tor_min_nom)/2; % [MPa]
tau_xz_tor_amp_nom =  (tau_xz_tor_max_nom - tau_xz_tor_min_nom)/2; % [MPa]


%%%%% Material data %%%%% Machine Design, Table A8 & A9 pg 1039-1040
material_key = [355, 4140];
material_data = struct('S_y', num2cell([355, 655]),...
                       'S_ut', num2cell([470, 758]));
material_table = dictionary(material_key, material_data);
S_y = material_table(material).S_y;
S_ut = material_table(material).S_ut;
S_yc = -S_y; % [MPa] Compressive yeild strength ! PLACEHOLDER VALUE must be incorporated with material data table


%%%%% Neubler's Constant for Steels %%%%% Machine Design, Table 6-6 page 382
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240]; % [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi_temp = S_ut * MPa_to_ksi; % [ksi] converts S_ut from MPa to ksi
S_ut_ksi = find_closest_value(S_ut_ksi_temp, S_ut_ksi_table);

a_sqrt_in = Neublers_table(S_ut_ksi); % [sqrt(in)]
a_sqrt_mm = a_sqrt_in * sqrt(25.4); % [sqrt(mm)] % MAS236, L3 s13 ~ Usikker på om dette er rett - RKH

% Notch sensitivity factor % MAS236 L3 s13 & Machine Design eq 6.13 pg 381
q = 1/(1 + (a_sqrt_mm/sqrt(r_fillet))); 


%%%%% Stress concentration factors %%%%% MAS236 L3 s12
% Geometrical (theoretical) stress concentration factors:
    % Tables/values: Machine Design Appendix C (page 1048-1049)
D_d_bend_key = [6.00, 3.00, 2.00, 1.50, 1.20, 1.10, 1.07, 1.05, ...
                1.03, 1.02, 1.01];
D_d_bend_values = struct('A', num2cell([0.87868, 0.89334, 0.90879, ...
                        0.93836, 0.97098, 0.95120, 0.97527, 0.98137, ...
                        0.98061, 0.96048, 0.91938]),...
                        'b', num2cell([-0.33243, -0.30860, -0.28598, ...
                        -0.25759, -0.21796, -0.23757, -0.20958, ...
                        -0.19653, -0.18381, -0.17711, -0.17032]));
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
D_d_axial_values = struct('A', num2cell([1.01470, 0.99957, 0.99682, ...
                            0.96272, 0.98084, 0.98450, 0.98498, ...
                            1.00480, 1.01220, 0.98413]), ...
                          'b', num2cell([-0.30035, -0.28221, -0.25751, ...
                          -0.25527, -0.22485, -0.20818, -0.19548, ...
                          -0.17076, -0.12474, -0.10474]));
D_d_axial_table = dictionary(D_d_axial_key, D_d_axial_values);

A_axial = D_d_axial_table(D_d).A;
b_axial = D_d_axial_table(D_d).b;

K_t_bend = A_bend*(r_fillet/d_shaft)^b_bend;    % Normal Stress - fig C-2
K_t_tor = A_tor*(r_fillet/d_shaft)^b_tor;       % Shear  Stress - fig C-3
K_t_axial = A_axial*(r_fillet/d_shaft)^b_axial; % Normal Stress - fig C-1

% For end-milled keyseat % Machine Design, fig 10-16, pg 615
if keyseat == 1
    K_t_tor = K_t_keyseat;
end

% Fatigue (dynamic) stress concentration foactors:
K_f_bend = 1 + q * (K_t_bend - 1);   % Normal Stress
K_f_tor = 1 + q * (K_t_tor - 1);     % Shear  Stress
K_f_axial = 1 + q * (K_t_axial - 1); % Normal Stress


%%%%% Mean & Amplitude stress w/stress concentration (Ductile materials) %%%%%
% MAS236 L3 s19 & s21

% Axial
sigma_x_axial_mean = sigma_x_axial_mean_nom * K_f_axial; % [MPa]
sigma_x_axial_amp =  sigma_x_axial_amp_nom  * K_f_axial; % [MPa]

% Shear
% tau_xy_shear_max = tau_xy_shear_max_nom * K_f_tor; % [MPa] Only xy !!
% tau_xy_shear_mean = tau_xy_shear_mean_nom * K_f_tor; % [MPa]
% tau_xy_shear_amp =  tau_xy_shear_amp_nom  * K_f_tor; % [MPa]
% tau_xz_shear_mean = tau_xz_shear_mean_nom * K_f_tor; % [MPa]
% tau_xz_shear_amp =  tau_xz_shear_amp_nom  * K_f_tor; % [MPa]

% Bending
% sigma_x_bend = sigma_x_bend_max_nom * K_f_bend; % [MPa]
sigma_x_bend_mean = sigma_x_bend_mean_nom * K_f_bend; % [MPa]
sigma_x_bend_amp =  sigma_x_bend_amp_nom  * K_f_bend; % [MPa]

% Torsion
% tau_xy_tor_max = tau_xy_tor_max_nom * K_f_tor; % [MPa]
tau_xy_tor_mean = tau_xy_tor_mean_nom * K_f_tor; % [MPa]
tau_xy_tor_amp =  tau_xy_tor_amp_nom  * K_f_tor; % [MPa]
tau_xz_tor_mean = tau_xz_tor_mean_nom * K_f_tor; % [MPa]
tau_xz_tor_amp =  tau_xz_tor_amp_nom  * K_f_tor; % [MPa]

% Resultant mean and amplitude % MAS236 L3 s23-24
sigma_x_mean = sigma_x_axial_mean + sigma_x_bend_mean; % [MPa]
sigma_x_amp =  sigma_x_axial_amp + sigma_x_bend_amp;   % [MPa]

% Shear neglected
tau_xy_mean = tau_xy_tor_mean; % [MPa]
tau_xy_amp =  tau_xy_tor_amp;  % [MPa]
tau_xz_mean = tau_xz_tor_mean; % [MPa]
tau_xz_amp =  tau_xz_tor_amp;  % [MPa]

% Included shear (done incorrectly)
% tau_xy_mean = tau_xy_shear_mean + tau_xy_tor_mean; % [MPa]
% tau_xy_amp =  tau_xy_shear_amp + tau_xy_tor_amp;   % [MPa]
% tau_xz_mean = tau_xz_shear_mean + tau_xz_tor_mean; % [MPa]
% tau_xz_amp =  tau_xz_shear_amp + tau_xz_tor_amp;   % [MPa]

% Von Mises Equivalent Stress
sigma_vm_mean = sqrt(sigma_x_mean^2 + 3*(tau_xy_mean^2+tau_xz_mean^2)); % [MPa]
sigma_vm_amp =  sqrt(sigma_x_amp^2 + 3*(tau_xy_mean^2+tau_xz_amp^2));   % [MPa]
sigma_vm_max = sigma_vm_mean + sigma_vm_amp; % [MPa]


%%%%% Correction factors %%%%%
% Load factor % MAS236 L3 s43 & Machine Design pg 366
C_load_table_key = ["Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial"];
C_load_table_value = [1 0.75 1 0.75 1];
C_load_table = dictionary(C_load_table_key, C_load_table_value);

C_load = C_load_table(load_type);

% Size factor % MAS236 L3 s44 & Machine Design pg 367
if d_shaft <= 8
    C_size = 1;
elseif d_shaft <= 250
    C_size = 1.189*d_shaft^(-0.097);
else
    C_size = 0.6;
end

% Surface factor % MAS236 L3 s45 & Machine Design pg 369
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

% Temperature factor % MAS236 L3 s47 & Machine Design pg 371
if operating_temperature <= 450
    C_temp = 1;
elseif operating_temperature <= 550
    C_temp = 1-0.0058*(operating_temperature-450);
else
    fprintf("Error: Operating temperature too high\n")
end

% Reliability factor % MAS236 L3 s48 & Machine Design pg 371
C_reliab_table_key = [50 90 95 99 99.9 99.99 99.999 99.9999];
C_reliab_table_value = [1 0.897 0.868 0.814 0.753 0.702 0.659 0.620];
C_reliab_table = dictionary(C_reliab_table_key, C_reliab_table_value);
C_reliab = C_reliab_table(reliability);

% Other factors: MAS236 lecture 3 s49

% For steels with "knee" % MAS236 L3 s39
if S_ut < 1400
    S_e_prime = 0.5*S_ut;
else
    S_e_prime = 700;
end


%%%%% Diameter Equation %%%%%
if first_iteration == "y"
    
    % Estimating stress geometric concentration factors 
                                % for preliminary stage % MAS236 L5 s10
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
    K_f_axial = K_t_axial;
    
    % Correction factors for preliminary stage % MAS236 L5 s12
    C_load = 1; % Bending
    C_size = 1;
end

% Endurance limit
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime; % MAS236 L4 s5

d_eq =  (((16*n_f)/pi) * (sqrt(4*(K_f_bend*M_amp)^2+3*(K_f_tor*T_amp)^2)/S_e + sqrt(4*(K_f_bend*M_mean)^2+3*(K_f_tor*T_mean)^2)/S_ut))^(1/3);
fprintf('d_eq = %.2f [mm],  Recomended shaft diameter\n', d_eq)

% Quick check: failure againt yels at the first cycle % MAS236 L5 s7
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

%%%%% Safety factors %%%%%
% Static safety factor
if sigma_vm_mean >= 0
    n_y = S_y/(sigma_vm_mean + sigma_vm_amp);
else
    n_y = S_y/abs(sigma_vm_mean - sigma_vm_amp);
end

if n_y > 1
    fprintf('n_y = %.2f --> No static failure\n', n_y)
else
    fprintf('n_y = %.2f --> Static failure\n', n_y)
end

% Fatigue safety factor
if sigma_vm_mean >= 0 && sigma_vm_mean + sigma_vm_amp < S_y
    sigma_rev = sigma_vm_amp/(1-(sigma_vm_mean/S_ut));
elseif sigma_vm_mean >= 0 && sigma_vm_mean + sigma_vm_amp > S_y
    sigma_rev = S_y;
elseif sigma_vm_mean < 0 && abs(sigma_vm_mean - sigma_vm_amp) < abs(S_yc)
    sigma_rev = sigma_vm_amp;
elseif sigma_vm_mean < 0 && abs(sigma_vm_mean - sigma_vm_amp) > abs(S_yc)
    sigma_rev = abs(S_yc);
else
    fprintf('\sigma_{rev} error in Modified Goodman Diagram')
end
n_f = S_e/sigma_rev;

if n_f > 1
    fprintf('n_f = %.2f --> No fatigue failure, infinete life\n', n_f)
else
    fprintf('n_f = %.2f --> Fatigue failure, finete life\n', n_f)
end


% Function to find the closest value rounded down in the list (Most conservative approach)
function closest_value = find_closest_value(value, list)
    list = list(list <= value); % Filter out values greater than the given value
    if isempty(list)
        error('No values in the list are less than or equal to the given value.');
    end
    [~, index] = min(abs(list - value)); % Find the closest value among the remaining ones
    closest_value = list(index);
end