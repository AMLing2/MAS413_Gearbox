% step : [mm] step to increase module at 
% sigma_b_lim : [MPa] material bending stress limit after safety factor 
% sigma_o_lim : [MPa] material conctact stress limit after safety factors 
% z : number of teeth on calculating gear
% n : [rpm] speed of gear 
% T : [Nmm] torque of gear 
% A : [m/s] operating factor, tab 2 pg 6 lec 4
% K_a : external dynamic factor, electric motor with light shock, tab 1 pg 5 lec 4
% lambda :  width factor, 8-12, pg 17 lec 1
% gamma : teeth form factor
%%%% sigma_o only variables:
% F_w : [sqrt(N/mm^2)] material factor
% F_c : edge form factor
% d_1 : pinion pitch diameter (use "pinion" if gear is pinion and d_1 is unknown)
% i : gear ratio for step
% Z_v : speed factor

function [m_n,sigma_b,sigma_o] = module_calc(increment, sigma_b_lim, ...
            sigma_o_lim, z, n, T, A, K_a, lambda, gamma, F_w, F_c, i)
    [m_n_b, sigma_b] = sigma_b_calc(increment, sigma_b_lim, z, n, T, A, K_a, lambda, gamma);
    [m_n_o, sigma_o] = sigma_o_calc(increment, sigma_o_lim, z, n, T, A, K_a, lambda, F_w, F_c, i);

    m_n = max(m_n_b, m_n_o);
end

function [m_n,sigma_b] = sigma_b_calc(increment, sigma_b_lim, z, n, T, A, K_a, lambda, gamma)
    m_n = increment; % start with step to prevent div by 0
    new_sigma_b = inf;
    while new_sigma_b > sigma_b_lim
        r = (m_n * z) / 2; % [mm] pitch circle radius
        V_t = n * ((2*pi)/60) * (r*1e-3); % [m/s] pitch speed
        K_V = (A + V_t) / A; % [-] inner dynamic factor, KGR lec 4 pg 6
        F_th = T/r; % [N] theoretical tangential force component, KGR lec 4 pg 5
        new_sigma_b = (F_th * K_a * K_V * gamma) / (m_n^2 * lambda); % [MPa], KGR lec 4 pg 5
    
        m_n = m_n + increment;
    end
    sigma_b = new_sigma_b;
end

function [m_n,sigma_o] = sigma_o_calc(increment, sigma_o_lim, z, n, T, A, ...
                                    K_a, lambda, F_w, F_c, i)
    m_n = increment; % start with step to prevent div by 0
    new_sigma_o = inf;
    Z_v_speed = [0.25,1,2,3,4,5,6,7,8,9,10,12]; % pg 10 lec 4 tab 4 [m/s]
    Z_v_factor = [0.835,0.842,0.855,0.877,0.905,0.932,0.960,0.980,1.0,...
        1.015,1.033,1.058];                     % pg 10 lec 4 tab 4 [-]
    Z_v_dic = dictionary(Z_v_speed,Z_v_factor);

    % initialize variables
    r = (m_n * z) / 2; % [mm]
    V_t = n * ((2*pi)/60) * (r * 1e-3); % pitch speed [m/s]
    Z_v = Z_v_dic( closest(Z_v_speed, V_t) ); % [-]

    while new_sigma_o > (sigma_o_lim*Z_v)
        d_1 = m_n*z; % [mm]
        r = (m_n * z) / 2; % [mm]
        V_t = n * ((2*pi)/60) * (r*1e-3); % pitch speed [m/s]
        K_V = (A + V_t) / A; % [-] dynamic factor
        F_th = T/r; % [N] theoretical tangential force component
        new_sigma_o = F_w * F_c * sqrt((F_th * K_a * K_V * (i+1))/ ...
                                     (m_n * lambda * d_1 * i )); % [MPa], KGR lec 4 pg 9

        Z_v = Z_v_dic(closest(Z_v_speed,V_t)); % from previous step
        m_n = m_n + increment;
    end
    sigma_o = new_sigma_o;
end