% chapter 10.12 machine design pg 620
%%%INPUTS:
% d_o : outer diameter of hub / gear [mm]
% d_s : outer diameter of shaft [mm]
% l : length of hub/gear engagement [mm]
% mu : coefficient of friction
% E_o : Young's modulus of gear / hub material [GPa]
% E_i : Young's modulus of shaft material [GPa]
% V_o : poisson's ratio of gear / hub material
% V_i : poisson's ratio of shaft material
% size_part : either "hub" or "shaft", changes relevant radius of chosen part

%%%OUTPUTS:
% p : pressure in the interference [Mpa]
% T_max : maximum torque the joint can withstand [Nm]
% r_c : changed radius of specified part (inner hole for hub, outer radius for shaft) [mm]
% h : gear hole diameter tolerence 
% s : shaft diameter tolerence 
% sigma_t_s = tangenial stress in shaft [Mpa]
% sigma_t_o = radial stress in shaft [Mpa]
% sigma_r_s = tangenial stress in hub / gear [Mpa]
% sigma_r_o = radial stress in hub / gear [Mpa]

function [p,T_max,r_c,h,s,sigma_t_s,sigma_t_o,sigma_r_s,sigma_r_o] = pressFits(d_o,d_s,l,mu ,E_o,E_i,V_o,V_i,size_part)
    r_i = 0; % [mm] solid shaft, no hollow inner diameter
    r_o = d_o/2; % [mm] hub / gear outer radius
    r = d_s/2; % [mm] shaft radius

    % class 8 fit (h8), Appendix E-1 fundamentals of machine component design pg 854
    C_h = 0.0052;
    C_s = 0.0052;
    C_i = 0.0010;
    h = C_h * d_s^(1/3); % gear hole diameter tolerence 
    s = C_s * d_s^(1/3); % shaft diameter tolerence 
    i = C_i * d_s; % [mm] average interference
    
    delta_r = i/2; % is dividing by 2 correct here?
    if size_part == "hub"
        r_c = r - delta_r; % [mm] radius of gear / hub inner hole
        if (r_c > r_o)
            warning("Shaft larger than gear")
        end
    elseif size_part == "shaft"
        r = r + delta_r;
        r_c = r; % [mm] radius of shaft
    else
        error('size_part must be either "shaft" or "hub"')
    end
    sigma = 2 * delta_r;

    p = 0.5*sigma / ...
        ( ((r/(E_o*1e3)) * ( ((r_o^2 + r^2)/(r_o^2 - r^2) ) + V_o)) + ...
          ((r/(E_i*1e3)) * ( ((r^2 + r_i^2)/(r^2 - r_i^2) ) - V_i))); % [Mpa] pressure in interface, eq 10.14a, pg 620 machine design
    T_max = (2*pi*r^2 * mu * p * l) * 1e-3; % [Nm] eq 10.14b pg 621 machine design

    % tangential and radial stresses for shaft and hub, eq 10.15a to 10.16b pg 621 machine design
    sigma_t_s = -p * ((r^2 + r_i^2) / (r^2 - r_i^2));
    sigma_t_o = -p * ((r_o^2 + r^2) / (r_o^2 - r^2));

    sigma_r_s = -p;
    sigma_r_o = -p;

    % stress concentrations: Figure 10-20 pg 622, no equation for it
    % l_d = l/d_s;
    % p_sigma = p / sigma_b;

end

