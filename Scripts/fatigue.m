% Conversion factors
MPa_to_ksi = 0.1450377377; % [MPa] to [ksi] conversion factor

%%%%% Imported values from shaftDesign.m %%%%%
% Forces
P = cs(1); % [N] % Axial
T = cs(2); % [Nmm] % Torque
M = cs(3); % [Nmm] % Bending
V_y = cs(4); % [N] % Shear
V_z = cs(5); % [N] % Shear

% Geometry for stress concentration
d_shaft  = cs(6); % [mm] Shaft diameter
%keyseat  = cs(7); % true -> keyseat, false -> shoulder-fillet
r_fillet = cs(8); % [mm] Notch fillet radius

if ~first_iteration
    % Calculated values
    R = (d_shaft/2);  % [mm] Shaft radius
    A = pi*R^2;       % [mm^2] Shaft area
    I = (pi/4)*R^4;   % [mm^4] Moment of inertia (I_x = I_y)
    I_p = (pi/2)*R^4; % [mm^4] Polar moment of inertia
end

%%%%% Material data %%%%%
% https://steelnavigator.ovako.com/steel-grades/s355j2/
S_y = S355J2(1); % Yield strength [MPa]
S_ut = S355J2(2); % Ultimate tensile sstrength [MPa]
E = S355J2(3); % Youngs module [Pa]
V_shaft = S355J2(4); % Poisson´s ratio


%%%%% Correction factors %%%%%
% Load factor % MAS236 L3 s43 & Machine Design pg 366
C_load_table_key = ["Pure bending" "Pure axial" "Pure torsion" ...
                    "Complex axial" "Complex non axial"]; % Key-Value Pair
C_load_table_value = [1 0.75 1 0.75 1];
C_load_table = dictionary(C_load_table_key, C_load_table_value);

C_load = C_load_table(load_type);

% Surface factor % MAS236 L3 s45 & Machine Design pg 369: for S_ut in [MPa]
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
    error("Operating temperature too high\n")
end

% Reliability factor % MAS236 L3 s48 & Machine Design pg 371
C_reliab_table_key = [50 90 95 99 99.9 99.99 99.999 99.9999]; % Key-Value
C_reliab_table_value = [1 0.897 0.868 0.814 0.753 0.702 0.659 0.620];
C_reliab_table = dictionary(C_reliab_table_key, C_reliab_table_value);
C_reliab = C_reliab_table(reliability);

% Other factors: MAS236 lecture 3 s49

% For steels with "knee" % MAS236 L3 s39
if S_ut < 1400 % [MPa]
    S_e_prime = 0.5*S_ut; % [MPa]
else
    S_e_prime = 700; % [MPa]
end


%%%%% Max & Min forces %%%%%
% Axial (constant)
P_max =  P; % [N]
P_min =  P; % [N]
P_m   = (P_max + P_min)/2; % [N]
P_a   = (P_max - P_min)/2; % [N]

% Shear (fully reversed)
V_y_max =  abs(V_y); % [N]
V_y_min = -V_y_max; % [N]
V_y_m   = (V_y_max + V_y_min)/2; % [N]
V_y_a   = (V_y_max - V_y_min)/2; % [N]
V_z_max =  abs(V_z); % [N]
V_z_min = -V_z_max; % [N]
V_z_m   = (V_y_max + V_y_min)/2; % [N]
V_z_a   = (V_y_max - V_y_min)/2; % [N]

% Bending (fully reversed)
M_max =  M; % [Nmm]
M_min = -M; % [Nmm]
M_m   = (M_max + M_min)/2; % [Nmm]
M_a   = (M_max - M_min)/2; % [Nmm]

% Torque (constant)
T_max =  T; % [Nmm]
T_min =  T; % [Nmm]
T_m   = (T_max + T_min)/2; % [Nmm]
T_a   = (T_max - T_min)/2; % [Nmm]


%%%%% Diameter Equation %%%%%
if first_iteration
    
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
    
    % Conservative estimate for preliminary stage (q is unknown)
    K_f_bend = K_t_bend;
    K_f_tor = K_t_tor;
    K_f_axial = K_t_axial;
    
    % Correction factors for preliminary stage % MAS236 L5 s12
    C_size = 1;

    % Endurance limit % MAS236 L4 s5
    S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime;
    
    % Diameter Equation % MAS236 L5 s6
    d_1st =  ( ( (16*n_f_desired) /pi) * ( sqrt(4*(K_f_bend*M_a)^2 + ...
                                       3*(K_f_tor*T_a)^2)/S_e + ...
                                  sqrt(4*(K_f_bend*M_m)^2 + ...
                                       3*(K_f_tor*T_m)^2)/S_ut) )^(1/3);
    fprintf('\n1st recomended shaft diameter --> d_1st = %.2f [mm]\n', d_1st)
    % Other formulation in equation (10.8) % Machine Design pg 600 & 653

    % Quick check: failure againt yield at the first cycle % MAS236 L5 s7
    % sigma_prime_amp = sqrt(((32*K_f_bend*M_amp)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_amp)/(pi*d_shaft^3)));
    % sigma_prime_mean = sqrt(((32*K_f_bend*M_mean)/(pi*d_shaft^3)) + 3*((16*K_f_tor*T_mean)/(pi*d_shaft^3))); 
    % sigma_max = sigma_prime_mean + sigma_prime_amp;
    % n_y = S_y / sigma_maxc
    return
end

%%%%% Neubler's Constant for Steels %%%%% Machine Design, Table 6-6 page 382
S_ut_ksi_table = [50 55 60 70 80 90 100 110 120 130 140 160 180 200 220 240]; % [ksi]
a_sqrt_in_table = [0.130 0.118 0.108 0.093 0.080 0.070 0.062 0.055 0.049 0.044 0.039 0.031 0.024 0.018 0.013 0.009]; % [in^1/2]
Neublers_table = dictionary(S_ut_ksi_table, a_sqrt_in_table);

S_ut_ksi_temp = S_ut * MPa_to_ksi; % [ksi] converts S_ut from MPa to ksi
S_ut_ksi = find_closest_value_down(S_ut_ksi_temp, S_ut_ksi_table);

a_sqrt_in = Neublers_table(S_ut_ksi); % [sqrt(in)]
a_sqrt_mm = a_sqrt_in * sqrt(25.4); % [sqrt(mm)] % MAS236, L3 s13

% Notch sensitivity factor % MAS236 L3 s13 & Machine Design eq 6.13 pg 381
q = 1/(1 + (a_sqrt_mm/sqrt(r_fillet)));


%%%%% Stress concentration factors %%%%% MAS236 L3 s12

% For end-milled keyseat % Machine Design, fig 10-16, pg 615
if cs(7)
    K_t_tor = cs(9); % Stress concentration for keyseat
   
    %%%%% Conservative estimate from 1st iteration for redundancy
    K_t_bend = 2.7;
    K_t_axial = 3.0;
    %%%%% (Should not apply)

    if any(abs([P, M, V_y, V_z]) > 1e-8)
    warning(['Check stress concentration factors K_t_axial & K_t_bend \n'...
             'P = %2f\nM = %2f\nV_y = %2f\nV_z = %2f\n'], P, M, V_y, V_z);
    end
else
    % Geometrical (theoretical) stress concentration factors for shoulder-fillet:
        % Tables/values: Machine Design Appendix C (page 1048-1049)
    
    % "Key-value pair" not physical/mechanical key
    D_d_bend_key = [6.00, 3.00, 2.00, 1.50, 1.20, 1.10, 1.07, 1.05, ...
                    1.03, 1.02, 1.01];
    D_d_bend_values = struct('A', num2cell([0.87868, 0.89334, 0.90879, ...
                            0.93836, 0.97098, 0.95120, 0.97527, 0.98137, ...
                            0.98061, 0.96048, 0.91938]),...
                            'b', num2cell([-0.33243, -0.30860, -0.28598, ...
                            -0.25759, -0.21796, -0.23757, -0.20958, ...
                            -0.19653, -0.18381, -0.17711, -0.17032]));
    D_d_bend_table = dictionary(D_d_bend_key, D_d_bend_values);
    D_d_bend = find_closest_value_up(cs(9), D_d_bend_key);
    A_bend = D_d_bend_table(D_d_bend).A;
    b_bend = D_d_bend_table(D_d_bend).b;
    
    % "Key-value pair" not physical/mechanical key
    D_d_tor_key = [2.00, 1.33, 1.20, 1.09];
    D_d_tor_values = struct('A', num2cell([0.86331, 0.84897, 0.83425, 0.90337]),...
                            'b', num2cell([-0.23865, -0.23161, -0.21649, -0.12692]));
    D_d_tor_table = dictionary(D_d_tor_key, D_d_tor_values);
    D_d_tor = find_closest_value_up(cs(9), D_d_tor_key);
    A_tor = D_d_tor_table(D_d_tor).A;
    b_tor = D_d_tor_table(D_d_tor).b;
    
    % "Key-value pair" not physical/mechanical key
    D_d_axial_key = [2.00, 1.50, 1.30, 1.20, 1.15, 1.10, 1.07, 1.05, 1.02, 1.01];
    D_d_axial_values = struct('A', num2cell([1.01470, 0.99957, 0.99682, ...
                                0.96272, 0.98084, 0.98450, 0.98498, ...
                                1.00480, 1.01220, 0.98413]), ...
                              'b', num2cell([-0.30035, -0.28221, -0.25751, ...
                              -0.25527, -0.22485, -0.20818, -0.19548, ...
                              -0.17076, -0.12474, -0.10474]));
    D_d_axial_table = dictionary(D_d_axial_key, D_d_axial_values);
    D_d_axial = find_closest_value_up(cs(9), D_d_axial_key);
    A_axial = D_d_axial_table(D_d_axial).A;
    b_axial = D_d_axial_table(D_d_axial).b;
    
    K_t_axial = A_axial*(r_fillet/d_shaft)^b_axial; % Normal Stress - fig C-1
    K_t_bend  = A_bend *(r_fillet/d_shaft)^b_bend;  % Normal Stress - fig C-2
    K_t_tor   = A_tor  *(r_fillet/d_shaft)^b_tor;   % Shear  Stress - fig C-3
end

% Fatigue (dynamic) stress concentration foactors:
K_f_bend  = 1 + q * (K_t_bend - 1);  % Normal Stress
K_f_tor   = 1 + q * (K_t_tor - 1);   % Shear  Stress
K_f_axial = 1 + q * (K_t_axial - 1); % Normal Stress


%%%%% Mean & Amplitude nominal stress %%%%% MAS236 L3 s15-18
% Axial
sigma_x_axial_max_nom = (P_max)/A; % [MPa]
sigma_x_axial_min_nom = (P_min)/A; % [MPa]
sigma_x_axial_m_nom   = (sigma_x_axial_max_nom + sigma_x_axial_min_nom)/2; % [MPa]
sigma_x_axial_a_nom   = (sigma_x_axial_max_nom - sigma_x_axial_min_nom)/2; % [MPa]

% Shear
tau_y_shear_max_nom = (4/3)*(V_y_max/A); % [MPa]
tau_y_shear_min_nom = (4/3)*(V_y_min/A); % [MPa]
tau_y_shear_m_nom   = (tau_y_shear_max_nom + tau_y_shear_min_nom)/2; % [MPa]
tau_y_shear_a_nom   = (tau_y_shear_max_nom - tau_y_shear_min_nom)/2; % [MPa]
tau_z_shear_max_nom = (4/3)*(V_z_max/A); % [MPa]
tau_z_shear_min_nom = (4/3)*(V_z_min/A); % [MPa]
tau_z_shear_m_nom   = (tau_z_shear_max_nom + tau_z_shear_min_nom)/2; % [MPa]
tau_z_shear_a_nom   = (tau_z_shear_max_nom - tau_z_shear_min_nom)/2; % [MPa]

% Bending
sigma_x_bend_max_nom = (M_max*R)/I; % [MPa]
sigma_x_bend_min_nom = (M_min*R)/I; % [MPa]
sigma_x_bend_m_nom   = (sigma_x_bend_max_nom + sigma_x_bend_min_nom)/2; % [MPa]
sigma_x_bend_a_nom   = (sigma_x_bend_max_nom - sigma_x_bend_min_nom)/2; % [MPa]

% Torsion
tau_tor_max_nom = (T_max*R)/I_p; % [MPa]
tau_tor_min_nom = (T_min*R)/I_p; % [MPa]
tau_tor_m_nom   = (tau_tor_max_nom + tau_tor_min_nom)/2; % [MPa]
tau_tor_a_nom   = (tau_tor_max_nom - tau_tor_min_nom)/2; % [MPa]

% Print stresses to compare
fprintf('tau_y_shear_max_nom  = %2f\n', tau_y_shear_max_nom)
fprintf('tau_z_shear_max_nom  = %2f\n', tau_z_shear_max_nom)
fprintf('tau_tor_max_nom      = %2f\n', tau_tor_max_nom)
fprintf('sigma_x_bend_max_nom = %2f\n', sigma_x_bend_max_nom)

%%%%% Mean & Amplitude stress w/stress concentration (Ductile materials) %%%%%
% MAS236 L3 s19 & s21

% Axial
sigma_x_axial_m = sigma_x_axial_m_nom * K_f_axial; % [MPa]
sigma_x_axial_a = sigma_x_axial_a_nom * K_f_axial; % [MPa]

% Shear
tau_y_shear_m = tau_y_shear_m_nom * K_f_tor; % [MPa]
tau_y_shear_a = tau_y_shear_a_nom * K_f_tor; % [MPa]
tau_z_shear_m = tau_z_shear_m_nom * K_f_tor; % [MPa]
tau_z_shear_a = tau_z_shear_a_nom * K_f_tor; % [MPa]

% Bending
sigma_x_bend_m = sigma_x_bend_m_nom * K_f_bend; % [MPa]
sigma_x_bend_a = sigma_x_bend_a_nom * K_f_bend; % [MPa]

% Torsion
tau_tor_m = tau_tor_m_nom * K_f_tor; % [MPa]
tau_tor_a = tau_tor_a_nom * K_f_tor; % [MPa]

% Resultant mean and amplitude % MAS236 L3 s23-24
sigma_x_m = sigma_x_axial_m + sigma_x_bend_m; % [MPa]
sigma_x_a = sigma_x_axial_a + sigma_x_bend_a; % [MPa]

% Shear neglected
% Von Mises Equivalent Stress
sigma_e_m   = sqrt(sigma_x_m^2 + 3*(tau_tor_m^2)); % [MPa]
sigma_e_a   = sqrt(sigma_x_a^2 + 3*(tau_tor_a^2)); % [MPa]
sigma_e_max = sigma_e_m + sigma_e_a; % [MPa]
fprintf('\nNeglected shear -->  sigma_e = %2f\n', sigma_e_max)

% Shear included
% Von Mises Equivalent Stress
sigma_e_m   = sqrt(sigma_x_m^2 + 3*(tau_tor_m^2+tau_y_shear_m^2+tau_z_shear_m^2)); % [MPa]
sigma_e_a   = sqrt(sigma_x_a^2 + 3*(tau_tor_a^2+tau_y_shear_a^2+tau_z_shear_a^2)); % [MPa]
sigma_e_max = sigma_e_m + sigma_e_a; % [MPa]
fprintf('Included shear  -->  sigma_e = %2f\n', sigma_e_max)


% Size factor % MAS236 L3 s44 & Machine Design pg 367
if d_shaft <= 8 % [mm]
    C_size = 1;
elseif d_shaft <= 250 % [mm]
    C_size = 1.189*d_shaft^(-0.097);
else
    C_size = 0.6;
end

% Endurance limit % MAS236 L4 s5
S_e = C_load*C_size*C_surf*C_temp*C_reliab*S_e_prime;

% Diameter Equation % MAS236 L5 s6
d_rec =  ( ( (16*n_f_desired) /pi) * ( sqrt(4*(K_f_bend*M_a)^2 + ...
                                   3*(K_f_tor*T_a)^2)/S_e + ...
                              sqrt(4*(K_f_bend*M_m)^2 + ...
                                   3*(K_f_tor*T_m)^2)/S_ut) )^(1/3);
fprintf('\nRecomended shaft diameter --> d_rec = %.2f [mm]\n', d_rec)
% Other formulation in equation (10.8) % Machine Design pg 600 & 653


%%%%% Safety factors %%%%%
% Static safety factor: MAS236 L5 s7
if sigma_e_m >= 0
    n_y = S_y/(sigma_e_m + sigma_e_a);
else
    n_y = S_y/abs(sigma_e_m - sigma_e_a);
end

if n_y > 1
    fprintf('\nNo static failure  --> n_y = %.2f\n', n_y)
else
    fprintf('\nStatic failure     --> n_y = %.2f\n', n_y)
end

% Fatigue safety factor: MAS236 L4 s19
if ( sigma_e_m >= 0 ) && ...
        ( (sigma_e_m + sigma_e_a) < S_y )
    sigma_rev = sigma_e_a / ( 1 - (sigma_e_m/S_ut) );
elseif ( sigma_e_m >= 0 ) && ...
        ( (sigma_e_m + sigma_e_a) > S_y )
    sigma_rev = S_y;
elseif ( sigma_e_m < 0 ) && ...
        ( abs(sigma_e_m - sigma_e_a) < abs(S_yc) )
    sigma_rev = sigma_e_a;
elseif ( sigma_e_m < 0 ) && ...
        ( abs(sigma_e_m - sigma_e_a) > abs(S_yc) )
    sigma_rev = abs(S_yc);
else
    error('\sigma_{rev} error in Modified Goodman Diagram')
end
n_f = S_e/sigma_rev;

if n_f > 1
    fprintf('No fatigue failure --> n_f = %.2f\n', n_f)
else
    fprintf('Fatigue failure    --> n_f = %.2f\n', n_f)
end


% Function to find the closest value rounded down in the list (Most conservative approach)
function closest_value_down = find_closest_value_down(value, list)
    list = list(list <= value); % Filter out values greater than the given value
    if isempty(list)
        error('No values in the list are less than or equal to the given value.');
    end
    [~, index] = min(abs(list - value)); % Find the closest value among the remaining ones
    closest_value_down = list(index);
end

% Function to find the closest value rounded up in the list (Most conservative approach)
function closest_value_up = find_closest_value_up(value, list)
    list = list(list >= value); % Filter out values less than the given value
    if isempty(list)
        error('No values in the list are greater than or equal to the given value.');
    end
    [~, index] = min(abs(list - value)); % Find the closest value among the remaining ones
    closest_value_up = list(index);
end