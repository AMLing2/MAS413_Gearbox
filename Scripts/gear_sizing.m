clc;close all;clear;

%TODO: 
% contact stress, usually sigma_o > sigma_b - done
% change m_n to m_t for all calcs - done
% check i_tot is 1% of 17.3 - done
% contact ratio between 1 and 2 - checked
% material factor see line 73
% stages must obey tables 12-4 12-5 pg 737, add check z2,z4 < 1309 ?
% increase lambda to 14?

% module of elasticity and material standards (no price):
% https://www.michael-smith-engineers.co.uk/mse/uploads/resources/useful-info/General-Info/MATERIAL-GRADE-COMPARISON-TABLE-for-Web.pdf
% material (diameter) (grade) (related link)
% Fe 430: (726mm) (S275) https://matweb.com/search/DataSheet.aspx?MatGUID=c2ba59bb365942a7b6da46f1cee370b8  
% Fe 590: (629mm) (S355) https://matweb.com/search/DataSheet.aspx?MatGUID=1dc0414bd1ea4061a5dc09382c455e2a
% C 45 N: (622mm) (AISI/ASTM 1045) https://matweb.com/search/DataSheet.aspx?MatGUID=2ca9b42e83894e8a8a61385fd7da63ae
% C 60 N: (579mm) (SAE/ASTM 1060) https://matweb.com/search/DataSheet.aspx?MatGUID=0a471605c1324daa910855e54a21fab3
% 34 Cr 4 V: (510mm) (AISI 5132/ 1.7033) https://matweb.com/search/DataSheet.aspx?MatGUID=4877d405464f448a96786c8cbd00d3b5
% 42 CrMo 4 V: (487mm) (ASTM A322) https://matweb.com/search/DataSheet.aspx?MatGUID=38108bfd64c44b4c9c6a02af78d5b6c6
% 16 MnCr 5: (285mm) (AISI 5115 / 1.7131) https://matweb.com/search/DataSheet.aspx?MatGUID=2ab813ffa05d40329dffe0ee7f58b5de
% 15 CrNi 6: (277mm) https://matweb.com/search/DataSheet.aspx?MatGUID=9ab3bf332758468ab36010790bd94349 ?

%%% Chosen Parameters %%%
material = "16 MnCr 5";
lambda = 14; % width factor, processed:  8-16, pg 17 lec 1, 14 to increase axial contact ratio to 1.15
n_f = 2; % safety factor
step = 0.05; % module increase for interference fit
%%% Chosen Parameters %%%

% Given Parameters:
n_1 = 1450; % [RPM]
P_1 = 12.5e3; % [W] input power
i_tot_og = 17.3;
alpha = 20;  % [deg] pressure angle, phi/theta
beta = 15;   % [deg] helix angle, psi

% Contact Ratio Minimum Requirement
mp_min = 1.4;  % machine design pg 738 - Contact Ratio
mF_min = 1.15; % machine design pg 796 - Transverse Contact Ratio

% Constant Factors from Tables
V_b = 1.7;  % [-] bending safety factor, KGR lec 4 pg 8
V_o = 1.25; % [-] contact safety factor, KGR lec 4 pg 10
K_L = 1.0;  % [-] lubrication factor, KGR lec 4 pg 10

% Min & max from Rollof: Figure 15-38 KGR Lecture 2 Slide 12
i1_min = 4.5;
i1_max = 5;
delta = 1e-3; % Iteration resolution
z_1 = 18; % Machine Design p737 Table 12-4: pressure angle
% zSmallest & zClosest are lists with format: [z_1,z_2,z_3,z_4,i_tot]
[zSmallest, zClosest] = grat2stage(i1_min, i1_max, delta, z_1, i_tot_og);
z_2 = zClosest(2);
z_3 = zClosest(3);
z_4 = zClosest(4);
i_tot = zClosest(5);

% gear ratios of stages
i_s1 = z_2/z_1;
i_s2 = z_4/z_3;

% check if new i_tot is within 1% of requirement
if ( abs(i_tot - i_tot_og)/i_tot_og ) > 0.01
    error("i_tot is greater than 1% of requirement")
end

% Gear Speed [RPM]
n_2 = n_1 / i_s1;
n_3 = n_2;
n_4 = n_3/ i_s2;

% Gear Velocity [rad/sec]
omega_1 = n_1 * (2*pi) / 60;
% omega_2 = n_2 * (2*pi) / 60;
% omega_3 = n_3 * (2*pi) / 60;
% omega_4 = n_4 * (2*pi) / 60;

% Torques in each gear [Nmm]
T_1 = P_1/omega_1 * 1e3;
T_2 = T_1 * i_s1;
T_3 = T_2;
T_4 = T_3 * i_s2;

%%%% add helical module calcs

% Material limits, KGR lec 4 pg 11 - table 5
sigma_b_lim_mat_list = [160,210,220,250,300,310,410,410]; % [MPa]
sigma_o_lim_mat_list = [430,520,540,610,715,760,1600,1900]; % [MPa]
E_mat_list = [200 200 210 210 205 205 200 210]*1e3; % [MPa]
mat_names = ["Fe 430", "Fe 590", "C 45 N", "C 60 N",...
    "34 Cr 4 V", "42 CrMo 4 V", "16 MnCr 5", "15 CrNi 6"];
sigma_b_lim_mat = dictionary(mat_names, sigma_b_lim_mat_list);
sigma_o_lim_mat = dictionary(mat_names, sigma_o_lim_mat_list);
E_mat_dic = dictionary(mat_names, E_mat_list);
Sy_mat = 417; % [Mpa] yield strenght of 16 MnCr 5, turn into list in the future? https://shop.machinemfg.com/composition-properties-and-uses-of-sae-aisi-5115-alloy-steel/

sigma_b_lim = sigma_b_lim_mat(material) / V_b;
% Z_v is calculated in the module_calc function
sigma_o_lim = (sigma_o_lim_mat(material) / V_o) * K_L;

%%% Bending stress sigma_b from lectures

A = 5; % [m/s] operating factor, tab 2 pg 6 lec 4
K_a = 1.25; % [-] external dynamic factor, electric motor with light shock, tab 1 pg 5 lec 4
F_w = sqrt(0.35 * E_mat_dic(material) ); % [sqrt(N/mm^2)] material factor, lec 4 pg 9
F_c = 1.76; % [-] edge form factor for alpha = 20, lec 4 pg 9

% gammas for gears
z = [20 25 30 40 60 80 100];
gammas = [2.9 2.73 2.60 2.45 2.30 2.24 2.21];
gamma_1 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
gamma_2 = tableInterpolation(z_2, z(5:6), gammas(5:6)); % teeth form factor, 79 teeth, tab 3 pg 7 lec 4
gamma_3 = 2.9; % teeth form factor, 18 teeth, tab 3 pg 7 lec 4
gamma_4 = tableInterpolation(z_4, z(5:6), gammas(5:6)); % teeth form factor, 71 teeth, tab 3 pg 7 lec 4

% calculating normal modules for each gear
increment = 0.01; % for module iteration
% T_1 = T_1*1e-3;
% T_2 = T_2*1e-3;
% T_3 = T_3*1e-3;
% T_4 = T_4*1e-3;

% initial guess of diameters
% d1 = increment*z1;
% d3 = increment*z3;

[m_n_1, sigma_b_1, sigma_o_1] = module_calc(increment, sigma_b_lim, ...
                                sigma_o_lim, z_1, n_1, T_1, A, ...
                                K_a, lambda, gamma_1, F_w, F_c, i_s1);
% d2 = m_n_1 * z_1;
% [m_n_2, sigma_b_2, sigma_o_2] = module_calc(increment, sigma_b_lim, sigma_o_lim, z_2, n_2, T_2, A, ...
        % K_a, lambda, gamma_2, F_w, F_c, d2, i_s1, z_1);
[m_n_3, sigma_b_3, sigma_o_3] = module_calc(increment, sigma_b_lim, ...
                                sigma_o_lim, z_3, n_3, T_3, A, ...
                                K_a, lambda, gamma_3, F_w, F_c, i_s2);
% d4 = m_n_3 * z_3;
% [m_n_4, sigma_b_4, sigma_o_4] = module_calc(increment, sigma_b_lim, sigma_o_lim, z_4, n_4, T_4, A, ...
%         K_a, lambda, gamma_4, F_w, F_c, d4, i_s2, z_3);

% d1 = m_n_1*z_1;
% d2 = m_n_1*z_2;
% d3 = m_n_3*z_3;
% d4 = m_n_3*z_4;
m_n_2 = m_n_1;
m_n_4 = m_n_3;

format long
% modules = table(m_n_1, m_n_2, m_n_3, m_n_4)
modules = table(m_n_1, m_n_3)
% sigmaB = table(sigma_b_1, sigma_b_2, sigma_b_3, sigma_b_4)
% sigmaO = table(sigma_o_1, sigma_o_2, sigma_o_3, sigma_o_4)
sigmaB = table(sigma_b_1, sigma_b_3)
sigmaO = table(sigma_o_1, sigma_o_3)
format short

% converting to transverse (helical) module
mt_1 = m_n_1 / cosd(beta);
mt_2 = m_n_2 / cosd(beta);
mt_3 = m_n_3 / cosd(beta);
mt_4 = m_n_4 / cosd(beta);
mt_s1 = max([mt_1,mt_2]); % 3.1180 w/ 15 CrNi 6
mt_s2 = max([mt_3,mt_4]); % 4.5334 w/ 15 CrNi 6

%%%%%%%% sizing calcs for helical gears
fits_unchecked = true;
stress_s1_checked = false;
while fits_unchecked
% pitch circle diameters iteration 1 [mm]
d_g1 = mt_s1 * z_1;
d_g2 = mt_s1 * z_2;
d_g3 = mt_s2 * z_3;
d_g4 = mt_s2 * z_4;

% width of helical gears [mm]
b_s1 = mt_s1 * lambda;
b_s2 = mt_s2 * lambda;

% top (ht) and bottom (hf) heights [mm]
ht_1 = mt_s1;
ht_2 = mt_s1; % see figure 12-8 pg 734 machine element
ht_3 = mt_s2; 
ht_4 = mt_s2; 
hf_1 = 1.25 * mt_s1;
hf_2 = 1.25 * mt_s1; 
hf_3 = 1.25 * mt_s2; 
hf_4 = 1.25 * mt_s2; 

% % pitch [mm]
% p_s1 = pi*mt_s1;
% p_s2 = pi*mt_s2;
% 
% % Diameter [mm], Juvinall eq (15.2) page 626
% d_g1 = p_s1*z_1/pi;
% d_g2 = p_s1*z_2/pi;
% d_g3 = p_s2*z_3/pi;
% d_g4 = p_s2*z_4/pi;

% addedum circle [mm]
dt_g1 = d_g1 + 2 * ht_1;
dt_g2 = d_g2 + 2 * ht_2;
dt_g3 = d_g3 + 2 * ht_3;
dt_g4 = d_g4 + 2 * ht_4;

% dedendum circle [mm]
df_g1 = d_g1 - 2 * hf_1;
df_g2 = d_g2 - 2 * hf_2;
df_g3 = d_g3 - 2 * hf_3;
df_g4 = d_g4 - 2 * hf_4;

% helical module iteration 2 and press fits:
if true%exist("shaftDesign.mat","file")
    if false %~exist("d_s_111","var") % only load once
        load("shaftDesign.mat")
    end
    d_shaft_g1 = 30; % [mm] % replace with loaded variable
    d_shaft_g2 = 50;
    d_shaft_g3 = 50;
    d_shaft_g4 = 60;

    mu = 0.175; % pg 621 machine design, between 0.15 and 0.2 for shrink fit hubs
    %mu = 0.74; % for static dry, mild steel on mild steel, tab 7-1 pg 464 machine design
    E_mat = E_mat_dic(material); % [MPa]
    V_mat = 0.29; % where does this come from?...
    % calculate press fits of each gear
    if ~stress_s1_checked
        % g1
        [p_g1,T_max_g1,d_i_g1,h_tol_g1,s_tol_g1,heat_temp_hub_g1,cool_temp_shaft_g1 ...
            ,sigma_t_s_g1,sigma_t_o_g1,sigma_r_s_g1,sigma_r_o_g1] = ... % stresses
            pressFitsShaft(df_g1,d_shaft_g1,b_s1,mu,E_mat*1e-3,E_mat*1e-3,V_mat,V_mat);
        sigma_t_o_g1
        sigma_r_o_g1
        
        % g2
        [p_g2,T_max_g2,d_i_g2,h_tol_g2,s_tol_g2,heat_temp_hub_g2,cool_temp_shaft_g2 ...
            ,sigma_t_s_g2,sigma_t_o_g2,sigma_r_s_g2,sigma_r_o_g2] = ... % stresses
            pressFitsShaft(df_g2,d_shaft_g2,b_s1,mu,E_mat*1e-3,E_mat*1e-3,V_mat,V_mat);
        sigma_t_o_g2
        sigma_r_o_g2
    
        % check if stresses are not too large for each gear:
        if (abs(sigma_t_s_g1) > (Sy_mat / n_f) &&  abs(sigma_r_o_g1) > (Sy_mat / n_f)) || ...
           (abs(sigma_t_s_g2) > (Sy_mat / n_f) &&  abs(sigma_r_o_g2) > (Sy_mat / n_f))  
            % increase module of first stage
            mt_s1 = mt_s1 + step;
            continue % restart loop
        else
            stress_s1_checked = true; % don't recalculate in future loops
        end
    end
    % check for stage 2
    % g3
    [p_g3,T_max_g3,d_i_g3,h_tol_g3,s_tol_g3,heat_temp_hub_g3,cool_temp_shaft_g3 ...
        ,sigma_t_s_g3,sigma_t_o_g3,sigma_r_s_g3,sigma_r_o_g3] = ... % stresses
        pressFitsShaft(df_g3,d_shaft_g3,b_s2,mu,E_mat*1e-3,E_mat*1e-3,V_mat,V_mat);
    % g4
    [p_g4,T_max_g4,d_i_g4,h_tol_g4,s_tol_g4,heat_temp_hub_g4,cool_temp_shaft_g4 ...
        ,sigma_t_s_g4,sigma_t_o_g4,sigma_r_s_g4,sigma_r_o_g4] = ... % stresses
        pressFitsShaft(df_g4,d_shaft_g4,b_s2,mu,E_mat*1e-3,E_mat*1e-3,V_mat,V_mat);

    % check
    if (abs(sigma_t_s_g3) > (Sy_mat / n_f) &&  abs(sigma_r_o_g3) > (Sy_mat / n_f)) || ...
       (abs(sigma_t_s_g4) > (Sy_mat / n_f) &&  abs(sigma_r_o_g4) > (Sy_mat / n_f))  
        % increase module of second stage
        mt_s2 = mt_s2 + step;
        continue % restart loop
    end
    if T_max_g1 < T_1*1e-3
        error("Shaft 1 cant transmit torque")
    elseif T_max_g2 < T_2*1e-3
        error("Shaft 2 cant transmit torque")
    elseif T_max_g3 < T_3*1e-3
        error("Shaft 3 cant transmit torque")
    elseif T_max_g4 < T_4*1e-3
        error("Shaft 4 cant transmit torque")
    end
    fits_unchecked = false; % exit loop
else
    warning("Unknown shaft diameters, skipping shrink fits")
    fits_unchecked = false; % dont loop
end
end
max_torques = table(T_max_g1, T_max_g2,T_max_g3,T_max_g4)
modules = table(mt_s1, mt_s2)

% % gearbox total length of gears
% l_tot = (dt_g1 + d_g2/2 + d_g3/2 + dt_g4)/1e3; % [m]


% % tooth thickness [mm]
% sn_s1 = p_s1/2 - 0.05 * mt_s1;
% sn_s2 = p_s2/2 - 0.05 * mt_s2;

% % hatch width [mm]
% en_s1 = p_s1/2 + 0.05 * mt_s1;
% en_s2 = p_s2/2 + 0.05 * mt_s2;
% 
% % axial pitch [mm] 
% px_s1 = p_s1 / sind(beta);
% px_s2 = p_s2 / sind(beta);
% 
% % diameteral pitch [mm]
% dp_s1 = pi/p_s1;
% dp_s2 = pi/p_s2;

% rough sum of material volume for gears [mm^3] -> [m^3]
volume_g1 = pi * b_s1 * (d_g1/2)^2 * 1e-9;
volume_g2 = pi * b_s1 * (d_g2/2)^2 * 1e-9;
volume_g3 = pi * b_s2 * (d_g3/2)^2 * 1e-9;
volume_g4 = pi * b_s2 * (d_g4/2)^2 * 1e-9;
material_sum = volume_g1 + volume_g2 + volume_g3 + volume_g4;
density = 7.85*1e3; % [kg/m^3]
mass_g1 = volume_g1 * density;
mass_g2 = volume_g2 * density;
mass_g3 = volume_g3 * density;
mass_g4 = volume_g4 * density;

% fillet radius
pd_in_s1 = 25.4/mt_s1; % [1/in] machine design eq 12.4d pg735
r_f_1 = (0.3/pd_in_s1)*25.4; %[in] -> [mm] machine design tab 12-1 pg735
pd_in_s2 = 25.4/mt_s2; % [1/in] machine design eq 12.4d pg735
r_f_2 = (0.3/pd_in_s2)*25.4; %[in] -> [mm] machine design tab 12-1 pg735

% center distance
C_s1 = d_g1/2 + d_g2/2; % [mm]
C_s2 = d_g3/2 + d_g4/2; % [mm]

%%% contact ratio %%%
Z_s1 =  sqrt((d_g1/2 + ht_1)^2 - (d_g1/2 * cosd(alpha))^2) + ...
        sqrt((d_g2/2 + ht_2)^2 - (d_g2/2 * cosd(alpha))^2) ...
        - C_s1 * sind(alpha); % eq 12.2 machine design pg 730
P_b_s1 = (pi * d_g1/ z_1) * cosd(alpha); % eq 12.3b machine design pg 734
CR_s1 = Z_s1/P_b_s1; % eq 12.7a machine design pg 738

Z_s2 =  sqrt((d_g3/2 + ht_3)^2 - (d_g3/2 * cosd(alpha))^2) + ...
        sqrt((d_g4/2 + ht_4)^2 - (d_g4/2 * cosd(alpha))^2) ...
        - C_s2 * sind(alpha); % eq 12.2 machine design pg 730
P_b_s1 = (pi * d_g3/ z_3) * cosd(alpha); % eq 12.3b machine design pg 734
CR_s2 = Z_s2/P_b_s1; % eq 12.7a machine design pg 738
if or( (CR_s1 < mp_min), (CR_s2 < mp_min) )
    error("Contact ratio is less than minimum")
end
%%% contact ratio %%%

%%% contact ratio - alternate from Lecture 3 page 12 %%%
    % More readable table from:
    % https://brainly.in/question/3177379
% Variable shift for consistency with table, copy existing
a_x = C_s1;                  % Center Distance
m_t = mt_s1;                 % Transverse Module
% alpha_t = atand( tand(alpha)/cosd(beta) );
alpha_t = alpha;             % Transverse (Helical) Angle
d_a1 = dt_g1;                % Addendum Diameter
d_a2 = dt_g2;                % Addendum Diameter
d_b1 = d_g1 * cosd(alpha_t); % Base Diameter
d_b2 = d_g2 * cosd(alpha_t); % Base Diameter
% From L3p12 - Common for all in table
alpha_wt = acosd( (d_b1 + d_b2) / (2*a_x) );
% From L3p12 - Contact Ratio of Helical Pair
epsilon_a = ( sqrt( (d_a1/2)^2 - (d_b1/2)^2 ) + ...
              sqrt( (d_a2/2)^2 - (d_b2/2)^2 ) - ...
              a_x * sind(alpha_wt) ) / ...
            ( pi*m_t*cosd(alpha_t) );
%%% contact ratio - alternate from Lecture 3 page 12 %%%
contactRatioStep1 = table(CR_s1, epsilon_a)

%%% axial contact ratio %%%
m_f_s1 = (b_s1 * tand(beta)) / (pi * mt_s1); % eq 13.5 machine design pg 796
m_f_s2 = (b_s2 * tand(beta)) / (pi * mt_s2);
if or( (m_f_s1 < 1), (m_f_s2 < 1) )
    error("Axial contact ratio is under 1")
elseif or((m_f_s1 < mF_min),(m_f_s2 < mF_min))
    warning("Axial contact ratio is low")
end
%%% axial contact ratio %%%

% Display results
T = table([d_g1; d_g2; d_g3; d_g4], [b_s1; b_s1; b_s2; b_s2], [mass_g1; mass_g2; mass_g3; mass_g4], ...
    'VariableNames', {'Diameter [mm]', 'Width [mm]', 'Mass [kg]'}, ...
    'RowNames', {'gear 1', 'gear 2', 'gear 3', 'gear 4'});

disp(T)

save("gear_sizes.mat")


function m_t = module_from_fit()
    
end