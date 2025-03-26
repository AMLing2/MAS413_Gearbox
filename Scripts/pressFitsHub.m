% chapter 10.12 machine design pg 620
% case: known hub inner diameter, increase shaft outer diamter to fit the hub
% (for bearings)
%%%INPUTS:
% d_h_o : outer diameter of hub / gear [mm]
% d_s : outer diameter of hub [mm]
% l : length of hub/gear engagement [mm]
% mu : coefficient of friction, static
% E_o : Young's modulus of gear / hub material [GPa]
% E_i : Young's modulus of shaft material [GPa]
% V_o : poisson's ratio of gear / hub material
% V_i : poisson's ratio of shaft material
% d_h_i : inner diameter of hub / bearing [mm]

%%%OUTPUTS:
% p : pressure in the interference [Mpa]
% T_max : maximum torque the joint can withstand [Nm]
% d_c : changed radius of outer shaft radius [mm]
% h : gear hole diameter tolerence
% s : shaft diameter tolerence
% heat_temp_hub : temperature of hub / bearing needed for the fit [deg C]
% temp_shaft : temperature of shaft needed for the fit [deg C]
% sigma_t_s = tangenial stress in shaft [Mpa]
% sigma_t_o = radial stress in shaft [Mpa]
% sigma_r_s = tangenial stress in hub / bearing [Mpa]
% sigma_r_o = radial stress in hub / bearing [Mpa]

function [p,T_max,d_c,h,s,heat_temp_hub,temp_shaft,sigma_t_s,sigma_t_o,sigma_r_s,sigma_r_o] ...
    = pressFitsHub(d_h_o,d_s,l,mu ,E_o,E_i,V_o,V_i,d_h_i)

    r_s_i = 0; % [mm] solid shaft, no hollow inner diameter
    r_h_o = d_h_o/2; % [mm] hub / gear outer radius

    % class 8 fit (h8), Appendix E-1 fundamentals of machine component design pg 854
    C_h = 0.0052;
    C_s = 0.0052;
    C_i_8 = 0.0010;
    C_i_7 = 0.0005;
    i_7 = C_i_7 * d_h_i; % [mm] average interference h7
    i_8 = C_i_8 * d_h_i; % [mm] average interference h8
    
    delta_r = i_7/2; % is dividing by 2 correct here?
    r = (d_h_i/2) + delta_r; % [mm] new radius of shaft
    d_c = r*2; % [mm] changed radius of shaft
    h = C_h * d_h_i^(1/3); % gear hole diameter tolerence 
    s = C_s * d_c^(1/3); % shaft diameter tolerence 

    % temperature calculations: pg 577 university physics book
    temp_room = 22; % room temperature [deg C]
    temp_shaft = temp_room;
    L_0 = 2*pi*(d_h_i/2); % initial circumference [mm]
    fit_tol = s+h; % extra radius for easier fit [mm]
    L_heated = 2*pi*((d_c + fit_tol)/2); % finishing length for fit [mm]
    alpha = 1.2e-5; % Coefficient of liear expansion [1/deg C], tab 17.1 pg 578
    heat_temp_hub = ((L_heated - L_0) / (alpha * L_0)) + temp_room; % [deg C] eq 17.6 pg 576

    max_heating = 100 + temp_room; % bad for bearings to be above 125 deg C
    if heat_temp_hub > max_heating
        %warning("bearing d = %dmm has to be cooled for interference fit",d_s)
        temp_shaft = temp_room - (heat_temp_hub - max_heating); % cool down shaft instead
        heat_temp_hub = max_heating;
    end

    % stress calculations
    sigma = 2 * delta_r;
    p = 0.5*sigma / ...
        ( ((r/(E_o*1e3)) * ( ((r_h_o^2 + r^2)/(r_h_o^2 - r^2) ) + V_o)) + ...
          ((r/(E_i*1e3)) * ( ((r^2 + r_s_i^2)/(r^2 - r_s_i^2) ) - V_i))); % [Mpa] pressure in interface, eq 10.14a, pg 620 machine design
    T_max = (2*pi*r^2 * mu * p * l) * 1e-3; % [Nm] eq 10.14b pg 621 machine design

    % tangential and radial stresses for shaft and hub, eq 10.15a to 10.16b pg 621 machine design
    sigma_t_s = -p * ((r^2 + r_s_i^2) / (r^2 - r_s_i^2));
    sigma_t_o = -p * ((r_h_o^2 + r^2) / (r_h_o^2 - r^2));

    sigma_r_s = -p;
    sigma_r_o = -p;

     if (d_s > d_c)
        warning("New shaft diameter is less than original")
        sigma_t_s = inf; % probably removing
        sigma_t_o = inf;
        sigma_r_s = inf;
        sigma_r_o = inf;
    end

    % stress concentrations: Figure 10-20 pg 622, no equation for it
    % l_d = l/d_s;
    % p_sigma = p / sigma_b;

end