% chapter 10.12 machine design pg 620
% case: known shaft diamter, set hub inner diamter to fit the shaft
%%%INPUTS:
% d_h_o : outer diameter of hub / gear [mm]
% d_s : outer diameter of shaft [mm]
% l : length of hub/gear engagement [mm]
% mu : coefficient of friction
% E_o : Young's modulus of gear / hub material [GPa]
% E_i : Young's modulus of shaft material [GPa]
% V_o : poisson's ratio of gear / hub material
% V_i : poisson's ratio of shaft material

%%%OUTPUTS:
% p : pressure in the interference [MPa]
% T_max : maximum torque the joint can withstand [Nm]
% d_c : changed diameter of inner radius of hub [mm]
% h : hub inner hole diameter tolerence
% s : shaft diameter tolerence
% heat_temp_hub : hub temperature needed for the fit [deg C]
% temp_shaft : shaft temperature needed (room temp if standard) [deg C]
% sigma_t_s = tangenial stress in shaft [MPa]
% sigma_t_h = radial stress in shaft [MPa]
% sigma_r_s = tangenial stress in hub / gear [MPa]
% sigma_r_h = radial stress in hub / gear [MPa]
% min_r_o = Minimum hub outer diamter [mm]

function [p,T_max,d_h_i,h,s,heat_temp_hub,temp_shaft,sigma_t_s,sigma_t_h,sigma_r_s,sigma_r_h] ...
    = shrinkFitGear(d_h_o,d_s,l,mu,E_o,E_i,V_o,V_i)

    r_i = 0; % [mm] solid shaft, no hollow inner diameter
    r_o = d_h_o/2; % [mm] hub / gear outer radius
    r = d_s/2; % [mm] shaft radius

    % class 8 fit (h7s6), Appendix E-1 fundamentals of machine component design pg 854
    C_h = 0.0052;
    C_s = 0.0052;
    C_i_8 = 0.0010;
    C_i_7 = 0.0005;
    i_7 = C_i_7 * d_s; % [mm] average interference class 7, h7p6
    i_8 = C_i_8 * d_s; % [mm] average interference class 8, h7s6
    % h7s6 from table 4.1a (Steinschaden Lecture 6 slide 12 of
        % "078-Engineering-Drawings-Lecture-Linear-Fits-Tolerances.pdf")
    
    delta_r = i_8/2;
    % resize:a
    d_h_i = d_s - i_8; % [mm] diameter of gear / hub inner hole
    h = C_h * d_h_i^(1/3); % hub inner hole diameter tolerence 
    s = C_s * d_s^(1/3); % shaft diameter tolerence 

    % temperature calculations: pg 577 university physics book
    temp_room = 22; % room temperature [deg C]
    temp_shaft = temp_room;
    L_0 = 2*pi*(d_h_i/2); % initial length (radius of outer - radius of inner) [mm]
    fit_tol = s+h; % extra radius for easier fit [mm]
    L_heated = 2*pi*((d_h_i + i_8 + fit_tol)/2); % finishing length for fit [mm]
    alpha = 1.2e-5; % Coefficient of linear expansion [1/deg C], tab 17.1 pg 578
    heat_temp_hub = ((L_heated - L_0) / (alpha * L_0)) + temp_room; % [deg C] eq 17.6 pg 576

    max_heating = 120; % [deg C] bad for bearings to be above 125 deg C
    if heat_temp_hub > max_heating
        %warning("Shaft d = %dmm has to be cooled for interference fit",d_s)
        temp_shaft = temp_room - (heat_temp_hub - max_heating); % cool down shaft instead
        heat_temp_hub = max_heating;
    end

    % stress calculations
    delta = 2 * delta_r;
    p = (0.5 * delta)/ ((r/(E_o*1e3))*((r_o^2 + r^2)/(r_o^2 - r^2) + V_o) + ((r/(E_i*1e3)) * ((r^2)/(r^2) -V_i)));
    % p = (0.5*sigma) / ...
    %     ( ((r/(E_o*1e3)) * ( ((r_o^2 + r^2)/(r_o^2 - r^2) ) + V_o)) + ...
    %       ((r/(E_i*1e3)) * ( ((r^2 + r_i^2)/(r^2 - r_i^2) ) - V_i))); % [MPa] pressure in interface, eq 10.14a, pg 620 machine design
    T_max = (2*pi*r^2 * mu * p * l) * 1e-3; % [Nm] eq 10.14b pg 621 machine design

    % tangential and radial stresses for shaft and hub, eq 10.15a to 10.16b pg 621 machine design
    sigma_t_s = -p * ((r^2 + r_i^2) / (r^2 - r_i^2));
    sigma_t_h =  p * ((r_o^2 + r^2) / (r_o^2 - r^2));

    sigma_r_s = -p;
    sigma_r_h = -p;
    if (d_h_i > d_h_o)
        warning("Shaft larger than gear")
        sigma_t_s = inf; % probably removing
        sigma_t_h = inf;
        sigma_r_s = inf;
        sigma_r_h = inf;
    end

    % n_f = 1.5; % saftey factor for minimum hub size calc;
    % sigma_y = 417; % [MPa]
    % min_r_o =  r*((E_o*n_f*r - E_i*n_f*r + 500*E_i*E_o*sigma*sigma_y + E_i*V_o*n_f*r - E_o*V_i*n_f*r)/(E_i*n_f*r + E_o*n_f*r - 500*E_i*E_o*sigma*sigma_y + E_i*V_o*n_f*r - E_o*V_i*n_f*r))^(1/2);

    % stress concentrations: Figure 10-20 pg 622, no equation for it
    % l_d = l/d_s;
    % p_sigma = p / sigma_b;

end