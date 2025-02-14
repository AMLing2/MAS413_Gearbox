% step : [mm] step to increase module at 
% sigma_b_lim : [MPa] material bending stress limit after safety factor 
% sigma_o_lim : [MPa] material conctact stress limit after safety factors 
% z : number of teeth on calculating gear
% n : [rpm] speed of gear 
% T : [Nmm] torque of gear 
% A : [m/s] operating factor, tab 2 pg 6 lec 4
% K_a : external dynamic factor, electric motor with light shock, tab 1 pg 5 lec 4
% lambda :  width factor, 8-12, pg 17 lec 1
% gamma : teeth form factor, 18 teeth, tab 3 pg 7 lec 4
%%%% sigma_o only variables:
% F_w : [sqrt(N/mm^2)] material factor
% F_c : edge form factor
% d_1 : pinion pitch diameter (use "pinion" if gear is pinion and d_1 is unknown)
% i : gear ratio for step

function [m_n,sigma_b,sigma_o] = module_calc(step, sigma_b_lim, sigma_o_lim, z, n, T, A, ...
        K_a, lambda, gamma, F_w, F_c, d_1, i)
    [m_n_b, sigma_b] = sigma_b_calc(step, sigma_b_lim, z, n, T, A, K_a, lambda, gamma);
    [m_n_o, sigma_o] = sigma_o_calc(step, sigma_o_lim, z, n, T, A, K_a, lambda, F_w, F_c, d_1, i);

m_n = max([m_n_b,m_n_o]);
end

function [m_n,sigma_b] = sigma_b_calc(step, sigma_b_lim, z, n, T, A, K_a, lambda, gamma)
    m_n = step; %start with step to prevent div by 0
    new_sigma_b = inf;
    while new_sigma_b > sigma_b_lim
        r = (m_n * z) / 2; % [mm] pitch circle radius
        V_t = n * ((2*pi)/60) * (r*1e-3); %pitch speed [m/s]
        K_V = (A + V_t) / A; % dynamic factor
        F_th = T/r; % theoretical tangential force component [N]
        new_sigma_b = (F_th * K_a * K_V * gamma) / (m_n^2 * lambda);
    
        m_n = m_n + step;
    end
    sigma_b = new_sigma_b;
end

function [m_n,sigma_o] = sigma_o_calc(step, sigma_o_lim, z, n, T, A, K_a, lambda, F_w, F_c, d_1, i)
    m_n = step; %start with step to prevent div by 0
    new_sigma_o = inf;
    d_1_new = d_1;
    while new_sigma_o > sigma_o_lim
        if strcmp(d_1,"pinion")
            d_1_new = m_n * z;
        end
        r = (m_n * z) / 2;
        V_t = n * ((2*pi)/60) * (r*1e-3); %pitch speed [m/s]
        K_V = (A + V_t) / A; % dynamic factor
        F_th = T/r; % theoretical tangential force component [N]
        new_sigma_o = F_w * F_c * sqrt((F_th * K_a * K_V * (i+1))/ ...
                                     (m_n * lambda * d_1_new * i ));

        m_n = m_n + step;
    end
    sigma_o = new_sigma_o;
end