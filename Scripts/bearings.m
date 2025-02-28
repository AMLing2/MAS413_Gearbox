clc;close all;clear;

%TODO:
% find actual datasheet

lifetime = 10; % [years]
work_cycle = 10; % [hours/day]

% teeth #
z_1 = 18;
z_2 = 79;
z_3 = 18;
z_4 = 71;

% gear ratios of stages
i_s1 = z_2/z_1;
i_s2 = z_4/z_3;
i_tot = i_s1 * i_s2;

% speed of gears [rpm]
n_1 = 1450; % input
n_2 = n_1 / i_s1;
n_3 = n_2;
n_4 = n_3/ i_s2;

% number of cycles through lifetime:
cycles_lifetime_sh1 = lifetime * 365.25 * work_cycle*60 * n_1
cycles_lifetime_sh2 = lifetime * 365.25 * work_cycle*60 * n_2
cycles_lifetime_sh3 = lifetime * 365.25 * work_cycle*60 * n_4

% Loads from mechOfMaterials_shaftx.m scripts:
% - radial [N]:
F_r_sh1 = 20;
F_r_sh2 = 20;
F_r_sh3 = 20;
% - axial  [N]:
F_a_sh1 = 5;
F_a_sh2 = 5;
F_a_sh3 = 5;

% minimum shaft diameters
d_min_sh1 = 0;
d_min_sh2 = 0;
d_min_sh3 = 0;

K_R = 0.62; % Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design
taper_ang = 20;
% bearing value lists from KRW
% https://www.krw.de/en/products/product-database/?tx_cskrwproducts_krwproducts%5Baction%5D=filter&tx_cskrwproducts_krwproducts%5Bcontroller%5D=Product&tx_cskrwproducts_krwproducts%5BcurrentPage%5D=15&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5Bcategories%5D%5B0%5D=31&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5BmaxInternalDiameter%5D=200&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5BminInternalDiameter%5D=10&cHash=feb6871c5b5e6249330c3cd6fdd08930

d_list = 60:10:150; % [mm]
bearing_name = ["32912" "32914" "32916" "32918" "32920" "32922" "32924" "32928" "32930"];
c_dyn_list = [49.5 76.8 80.6 105 136 143 185 222 228 305] * 1e3;% [kN] -> [N]
c_st_list = [82.3 123 136 178 231 253 326 428 545] * 1e3; % [kN] -> [N]

[b_index_1,l1] = tapered_bearing_sizing(d_min_sh1 ,F_r_sh1,F_a_sh1,cycles_lifetime_sh1, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list)
[b_index_2,l2] = tapered_bearing_sizing(d_min_sh2 ,F_r_sh1,F_a_sh1,cycles_lifetime_sh1, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list);
[b_index_3,l3] = tapered_bearing_sizing(d_min_sh3 ,F_r_sh1,F_a_sh1,cycles_lifetime_sh1, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list);

function [bearing_index,lifetime] = tapered_bearing_sizing(d_min,F_r,F_a,cycles,taper_ang,K_R,d_list,C_dyn_list,C_st_list)
% single row bearing

% d_min: minimum rod diameter [mm] % add max diameter?
% F_r: Force in the radial direction [N]
% F_a: Force in the axial direction [N]
% cycles: minimum lifetime cycles
% taper_ang: taper angle of bearing, must be 20, 25, 30, 35, or 40 [deg]
% K_R: Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design
% d_list: list of bearing inner diameter [mm]
% C_dyn_list: list of bearing dynamic load rating [N]
% C_st_list: list of bearing static load rating [N]

    V = 1.0; % rotation factor, rotating inner ring
    if taper_ang == 20 % values from Fig 11-24 pg 705 machine design
        X = 0.43; % radial factor
        Y = 1; % thrust (axial) factor
        e = 0.57; % axial/thrust factor minimum
    elseif taper_ang == 25
        X = 0.41;
        Y = 0.87;
        e = 0.68;
    elseif taper_ang == 30
        X = 0.39;
        Y = 0.76;
        e = 0.80;
    elseif taper_ang == 35
        X = 0.37;
        Y = 0.66;
        e = 0.95;
    elseif taper_ang == 40
        X = 0.35;
        Y = 0.57;
        e = 1.14;
    else
        error("taper angle needs to be 20, 25, 30, 35, or 40 deg")
    end
    if (F_a/(F_r*V)) <= e
        X = 1; % eq. 11.22b
        Y = 0;
    end
    P = X * V * F_r + Y * F_a; % [N] equivalent load, eq 11.22a pg 704 machine design
    bearing_index = -1;
    lifetime = -1;
    for i = 1:length(d_list)
        if d_list(i) > d_min
            L_10 = K_R*(C_dyn_list(i)/P)^(10/3);
            if (L_10 > cycles) && (P < C_st_list(i))
                bearing_index = i;
                lifetime = L_10;
                break
            end
        end
    end
    if bearing_index == -1
        error("No suitable bearing found")
    end
end