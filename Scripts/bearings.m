clc; close all; clear;

lifetime = 10; % [years]
daysPerYear = 365.25; % [days/year]
work_cycle = 10; % [hours/day]
minPerHour = 60; % [min/hour]
ly = lifetime * daysPerYear * work_cycle * minPerHour; % [min]

% Import from Gear Sizing
load('gear_sizes.mat', 'z_1','z_2','z_3','z_4','i_tot','i_s1','i_s2', ...
    'n_1','n_2','n_3','n_4')

% Import from Loading Diagrams
load('loadingDiagram_shaft1.mat', 'xz_P', 'xy_V', 'xz_V')
shaft1_Fa = xz_P;
shaft1_Fr = sqrt(xy_V^2 + xz_V^2);
load('loadingDiagram_shaft2.mat', 'xz_P', 'xy_V', 'xz_V')
shaft2_Fa = xz_P;
shaft2_Fr = sqrt(xy_V^2 + xz_V^2);
load('loadingDiagram_shaft3.mat', 'xz_P', 'xy_V', 'xz_V')
shaft3_Fa = xz_P;
shaft3_Fr = sqrt(xy_V^2 + xz_V^2);

% Import from Shaft Design
load('shaftDesign.mat', 'd_S11', 'd_C')

% number of cycles through lifetime:
cycles_lifetime_sh1 = ly * n_1
cycles_lifetime_sh2 = ly * n_2
cycles_lifetime_sh3 = ly * n_4

% Reliability factor, weibull distribution - tab 11-5 pg 701 machine design
K_R = 0.62; % R% = 95

% rough calculation of angle since not specified in datasheets
op = (200-150)/4; 
h = sqrt(op^2 + 44^2);
theta = asind(op/h);
taper_ang_list = 20:5:40;
[~,ang_index] = min(taper_ang_list-theta);
taper_ang = taper_ang_list(ang_index);
% bearing value lists from KRW
% https://www.krw.de/en/products/product-database/?tx_cskrwproducts_krwproducts%5Baction%5D=filter&tx_cskrwproducts_krwproducts%5Bcontroller%5D=Product&tx_cskrwproducts_krwproducts%5BcurrentPage%5D=15&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5Bcategories%5D%5B0%5D=31&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5BmaxInternalDiameter%5D=200&tx_cskrwproducts_krwproducts%5BfilterRequest%5D%5BminInternalDiameter%5D=10&cHash=feb6871c5b5e6249330c3cd6fdd08930

d_list = 50:10:150; % [mm]
bearing_name = ["32912" "32914" "32916" "32918" "32920" "32922" "32924" "32928" "32930"];
c_dyn_list = [49.5 76.8 80.6 105 136 143 185 222 228 305] * 1e3; % [N]
c_st_list = [82.3 123 136 178 231 253 326 428 545] * 1e3; % [N]

[b_index_1,n1] = tapered_bearing_sizing(d_min_sh1 ,F_r_sh1,F_a_sh1,cycles_lifetime_sh1, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list)
[b_index_2,n2] = tapered_bearing_sizing(d_min_sh2 ,F_r_sh1,F_a_sh2,cycles_lifetime_sh2, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list);
[b_index_3,n3] = tapered_bearing_sizing(d_min_sh3 ,F_r_sh1,F_a_sh3,cycles_lifetime_sh3, ...
    taper_ang, K_R, d_list, c_dyn_list, c_st_list);

% Deep Groove Ball Bearing Selection
% need to update c_dyn_list and c_st_list with vales from:
% https://www.skf.com/group/products/rolling-bearings/ball-bearings/deep-groove-ball-bearings#cid-493604
d_list_ball = 50:5:150; % [mm]
[b_index_m,n_bm] = ball_bearing_sizing(d_min_sh1_m ,F_r_sh1,0,cycles_lifetime_sh1, ...
    K_R, d_list_ball, c_dyn_list, c_st_list); 
[b_index_p,n_bp] = ball_bearing_sizing(d_min_sh3_p ,F_r_sh3,0,cycles_lifetime_sh3, ...
    K_R, d_list_ball, c_dyn_list, c_st_list); 


%% Tapered Roller Bearing Selection
function [bearing_index,lifetime] = tapered_bearing_sizing(d_min,F_r,F_a,cycles,taper_ang,K_R,d_list,C_dyn_list,C_st_list)
    % Single Row Tapered Roller Bearing
% d_min: minimum rod diameter [mm] % add max diameter?
% F_r: Force in the radial direction [N]
% F_a: Force in the axial direction [N]
% cycles: minimum lifetime cycles
% K_R: Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design
% d_list: list of bearing inner diameter [mm]
% C_dyn_list: list of bearing dynamic load rating [N]
% C_st_list: list of bearing static load rating [N]

    V = 1.0; % rotation factor, rotating inner ring
    if taper_ang == 20 % values from Fig 11-24 pg 705 machine design
        X = 0.43; % radial factor
        Y = 1; % thrust (axial) factor
        e = 0.57; % limiting factor
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
        %fprintf("ignoring axial factor\n")
        X = 1; % eq. 11.22b
        Y = 0;
    end
    P = X * V * F_r + Y * F_a; % [N] equivalent load, eq 11.22a pg 704 machine design
    bearing_index = -1;
    lifetime = -1;
    for i = 1:length(d_list)
        if d_list(i) > d_min
            L_P = ( K_R*(C_dyn_list(i)/P)^(10/3) ) * 1e6; % eq 11.20b pg 701 machine design
            if (L_P > cycles) && (P < C_st_list(i)) % is (P < C_st_list(i)) valid?
                bearing_index = i;
                lifetime = L_P;
                break
            end
        end
    end
    if bearing_index == -1
        error("No suitable bearing found")
    end
end


%% Ball Bearing Selection
function [bearing_index,lifetime] = ball_bearing_sizing(d_min,F_r,F_a,cycles,K_R,d_list,C_dyn_list,C_st_list)
    % Single Row Deep Groove (Conrad) Ball Bearing
% d_min: minimum rod diameter [mm] % add max diameter?
% F_r: Force in the radial direction [N]
% F_a: Force in the axial direction [N]
% cycles: minimum lifetime cycles
% K_R: Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design
% d_list: list of bearing inner diameter [mm]
% C_dyn_list: list of bearing dynamic load rating [N]
% C_st_list: list of bearing static load rating [N]

    % values from Fig 11-24 pg 705 machine design:
    V = 1.0; % rotation factor, rotating inner ring
    X = 0.56; % radial factor
    Y_list = [2.3,1.99,1.71,1.55,1.45,1.31,1.15,1.04,1.00];
    Fa_C0_list = [0.014,0.028,0.056,0.084,0.11,0.17,0.28,0.42,0.56];
    e_list = [0.19, 0.22, 0.26, 0.28, 0.30, 0.34, 0.38, 0.42, 0.44];
    
    bearing_index = -1;
    lifetime = -1;
    for i = 1:length(d_list)
        if d_list(i) > d_min
            fa_c0 = F_a / C_st_list(i);
            [~,Y_index] = closest(Fa_C0_list,fa_c0);
            Y = Y_list(Y_index);
            e = e_list(Y_index);
    
            if (F_a/(F_r*V)) <= e
                %fprintf("ignoring axial factor\n")
                X = 1; % eq. 11.22b
                Y = 0;
            end
            P = X * V * F_r + Y * F_a; % [N] equivalent load, eq 11.22a pg 704 machine design

            L_P = ( K_R*(C_dyn_list(i)/P)^3 ) * 1e6; % eq 11.20a pg 701 machine design
            if (L_P > cycles) && (P < C_st_list(i))
                bearing_index = i;
                lifetime = L_P;
                break
            end
        end
    end
    if bearing_index == -1
        error("No suitable bearing found")
    end
end

%{
 Find corresponding index:
    https://www.mathworks.com/matlabcentral/answers/
    152301-find-closest-value-in-array#answer_1559017
%}
function [cl,closestIndex] = closest(arr,val) 
    [~,closestIndex] = min(arr-val.', [], ComparisonMethod = "abs");
    cl = arr(closestIndex);
end