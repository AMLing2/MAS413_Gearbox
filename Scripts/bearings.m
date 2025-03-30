clc; close all; clear;

lifetime = 5; % [years]
daysPerYear = 260; % [days/year] work days
work_cycle = 10; % [hours/day]
minPerHour = 60; % [min/hour]
hours = lifetime * daysPerYear * work_cycle % [hours]
ly = lifetime * daysPerYear * work_cycle * minPerHour; % [min]

% Import shaft speeds from Gear Sizing
if exist(fullfile("export_import","gear_sizes.mat"),'file')
    load(fullfile("export_import","gear_sizes.mat"),"i_tot","i_s1","i_s2", ...
        "n_1","n_2","n_3","n_4") % [rpm]
else
    error("Run gear_sizing.m first")
end

% Import from Loading Diagrams: Axial Load Fa & Radial Load Fr
if exist(fullfile("export_import","loadingDiagram_shaft1.mat"),'file')
    load(fullfile("export_import","loadingDiagram_shaft1.mat"), "B_Fa", "B_Fr", "C_Fa", "C_Fr")
    load(fullfile("export_import","loadingDiagram_shaft2.mat"), "D_Fa", "D_Fr", "E_Fa", "E_Fr")
    load(fullfile("export_import","loadingDiagram_shaft3.mat"), "F_Fa", "F_Fr", "G_Fa", "G_Fr")
else
    error("Run loadingDiagrams.m first")
end

% Import diameters from Shaft Design
if exist(fullfile("export_import","shaft_design.mat"),'file') 
    load(fullfile("export_import","shaft_design.mat"), "aaaa")
else % initial run
    warning("Run shaftDesign.m first") % change to error in the future
    % First run initial diamterers
    d_min_B = 24; % [mm]
    d_min_C = 24; % [mm]
    d_min_D = 0; % [mm]
    d_min_E = 0; % [mm]
    d_min_F = 65; % [mm]
    d_min_G = 66; % [mm]
end

% load bearing .CSV file
% data from: https://www.skf.com/group/products/rolling-bearings/ball-bearings/deep-groove-ball-bearings#cid-493604
% csvfile = "../Data/combined_ballBearings_manual2.csv";
csvfile = fullfile("..","Data","combined_ballBearings_manual2.csv");
b_data = readtable(csvfile,"NumHeaderLines",9,"DecimalSeparator",".","Delimiter",";");
b_data.Properties.VariableNames = ["num","d","D","B","C","C0","Pu",...
    "ref_speed","max_speed","mass","name_null","designations", ...
    "capped_one_side","D_null","d1","d2","D1","D2","r1,2","da_min", "da_max","Da_max","ra_max","kr","f0"];

% number of cycles through lifetime:
cycles_lifetime_sh1 = ly * n_1; % number of cycles
cycles_lifetime_sh2 = ly * n_2;
cycles_lifetime_sh3 = ly * n_4;
cycles_tab = table(cycles_lifetime_sh1,cycles_lifetime_sh2,cycles_lifetime_sh3)

% Reliability factor, weibull distribution - tab 11-5 pg 701 machine design
K_R = 1.0; % R% = 90

% shaft 1 bearings B , C
[b_index_B,n_bB] = ball_bearing_sizing(d_min_B ,B_Fr,B_Fa,cycles_lifetime_sh1, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side); 
[b_index_C,n_bC] = ball_bearing_sizing(d_min_C ,C_Fr,C_Fa,cycles_lifetime_sh1, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side); 
% shaft 2: D E
[b_index_D,n_bD] = ball_bearing_sizing(d_min_D ,D_Fr,D_Fa,cycles_lifetime_sh2, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side); 
[b_index_E,n_bE] = ball_bearing_sizing(d_min_E ,E_Fr,E_Fa,cycles_lifetime_sh2, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side); 
% shaft 3 : F G
[b_index_G,n_bG] = ball_bearing_sizing(d_min_G ,G_Fr,G_Fa,cycles_lifetime_sh3, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side);
[b_index_F,n_bF] = ball_bearing_sizing(d_min_F ,F_Fr,F_Fa,cycles_lifetime_sh3, ...
    K_R, b_data.d, b_data.C, b_data.C0,b_data.f0,b_data.capped_one_side);

% display as tables:

data_names = ["Bore diameter [mm]";
    "Outer diameter [mm]";
    "Width [mm]";
    "da_min [mm]";
    "da_max [mm]";
    "ra_max [mm]"
    "mass [kg]";
    "limiting speed [rpm]";
    "Name";
    "lifetime cycles [millions of revs]";
    "lifetime hours [h]"];

% shaft 1 B C
% bearing B
lifetime_hour_B = (n_bB / n_1) / minPerHour;
b_B = b_data.B(b_index_B); % width of bearing
fillet_r_B = b_data.ra_max(b_index_B); % save maximum shaft fillet radius [mm]
data_bearing_B = [b_data.d(b_index_B);
    b_data.D(b_index_B);
    b_B;
    b_data.da_min(b_index_B);
    b_data.da_max(b_index_B);
    b_data.ra_max(b_index_B);
    b_data.mass(b_index_B);
    b_data.max_speed(b_index_B);
    string(b_data.capped_one_side(b_index_B));
    n_bB*1e-6;
    lifetime_hour_B];
bearing_tab_B = table(data_bearing_B,data_names,'VariableNames',["Bearing B data","variables"]');
disp(bearing_tab_B);

% bearing C
lifetime_hour_C = (n_bC / n_1) / minPerHour;
b_C = b_data.B(b_index_C);
fillet_r_C = b_data.ra_max(b_index_C);
data_bearing_C = [b_data.d(b_index_C);
    b_data.D(b_index_C);
    b_C;
    b_data.da_min(b_index_C);
    b_data.da_max(b_index_C);
    b_data.ra_max(b_index_C);
    b_data.mass(b_index_C);
    b_data.max_speed(b_index_C);
    string(b_data.capped_one_side(b_index_C));
    n_bC*1e-6;
    lifetime_hour_C];
bearing_tab_C = table(data_bearing_C,data_names,'VariableNames',["Bearing C data","variables"]');
disp(bearing_tab_C);

% shaft 2 bearings D E

% bearing D
lifetime_hour_D = (n_bD / n_2) / minPerHour;
b_D = b_data.B(b_index_D);
fillet_r_D = b_data.ra_max(b_index_D);
data_bearing_D = [b_data.d(b_index_D);
    b_data.D(b_index_D);
    b_data.B(b_index_D);
    b_data.da_min(b_index_D);
    b_data.da_max(b_index_D);
    b_data.ra_max(b_index_D);
    b_data.mass(b_index_D);
    b_data.max_speed(b_index_D);
    string(b_data.capped_one_side(b_index_D));
    n_bD*1e-6;
    lifetime_hour_D];
bearing_tab_D = table(data_bearing_D,data_names,'VariableNames',["Bearing D data","variables"]');
disp(bearing_tab_D);

% bearing E
lifetime_hour_E = (n_bE / n_2) / minPerHour;
b_E = b_data.B(b_index_E);
fillet_r_E = b_data.ra_max(b_index_E);
data_bearing_E = [b_data.d(b_index_E);
    b_data.D(b_index_E);
    b_data.B(b_index_E);
    b_data.da_min(b_index_E);
    b_data.da_max(b_index_E);
    b_data.ra_max(b_index_E);
    b_data.mass(b_index_E);
    b_data.max_speed(b_index_E);
    string(b_data.capped_one_side(b_index_E));
    n_bE*1e-6;
    lifetime_hour_E];
bearing_tab_E = table(data_bearing_E,data_names,'VariableNames',["Bearing E data","variables"]');
disp(bearing_tab_E);

% shaft 3 F G
% bearing F
lifetime_hour_F = (n_bF / n_1) / minPerHour;
b_F = b_data.B(b_index_F);
fillet_r_F = b_data.ra_max(b_index_F);
data_bearing_F = [b_data.d(b_index_F);
    b_data.D(b_index_F);
    b_data.B(b_index_F);
    b_data.da_min(b_index_F);
    b_data.da_max(b_index_F);
    b_data.ra_max(b_index_F);
    b_data.mass(b_index_F);
    b_data.max_speed(b_index_F);
    string(b_data.capped_one_side(b_index_F));
    n_bF*1e-6;
    lifetime_hour_F];
bearing_tab_F = table(data_bearing_F,data_names,'VariableNames',["Bearing F data","variables"]');
disp(bearing_tab_F);

% bearing G
lifetime_hour_G = (n_bG / n_1) / minPerHour;
b_G = b_data.B(b_index_G);
fillet_r_G = b_data.ra_max(b_index_G);
data_bearing_G = [b_data.d(b_index_G);
    b_data.D(b_index_G);
    b_data.B(b_index_G);
    b_data.da_min(b_index_G);
    b_data.da_max(b_index_G);
    b_data.ra_max(b_index_G);
    b_data.mass(b_index_G);
    b_data.max_speed(b_index_G);
    string(b_data.capped_one_side(b_index_G));
    n_bG*1e-6;
    lifetime_hour_G];
bearing_tab_G = table(data_bearing_G,data_names,'VariableNames',["Bearing G data","variables"]');
disp(bearing_tab_G);

% Calculate shrink fits for each bearing
% B C
[d_B,h_B,s_B,temp_B,temp_shaft_B] = ...
    shrinkFitBearing(b_data.D(b_index_B),b_data.d(b_index_B));
[d_C,h_C,s_C,temp_C,temp_shaft_C] = ...
    shrinkFitBearing(b_data.D(b_index_C),b_data.d(b_index_C));
% E D
[d_E,h_E,s_E,temp_E,temp_shaft_E] = ...
    shrinkFitBearing(b_data.D(b_index_E),b_data.d(b_index_E));
[d_D,h_D,s_D,temp_D,temp_shaft_D] = ...
    shrinkFitBearing(b_data.D(b_index_D),b_data.d(b_index_D));
% F G
[d_F,h_F,s_F,temp_F,temp_shaft_F] = ...
    shrinkFitBearing(b_data.D(b_index_F),b_data.d(b_index_F));
[d_G,h_G,s_G,temp_G,temp_shaft_G] = ...
    shrinkFitBearing(b_data.D(b_index_G),b_data.d(b_index_G));

% print shaft diameters
table(d_C,d_B,d_E,d_G,d_F)

% saving data:
clear b_data % remove table before saving
save(fullfile("export_import","bearings.mat"))

% save only diameters
save(fullfile("export_import","diameter_shaft_bearings.mat"),"d_B","d_C", ...
    "d_E", "d_D", "d_F", "d_G");

%% Ball Bearing Selection
function [bearing_index,lifetime] = ball_bearing_sizing(d_min,F_r,F_a,cycles,K_R,d_list,C_dyn_list,C_0_list,f0_list,desg_list)
    % Single Row Deep Groove (Conrad) Ball Bearing
% d_min: minimum rod diameter [mm] % add max diameter?
% F_r: Force in the radial direction [N]
% F_a: Force in the axial direction [N]
% cycles: minimum lifetime cycles [cycles]
% K_R: Reliability factor for weibull distribution, tab 11-5 pg 7-1 machine design
% d_list: list of bearing inner diameter [mm]
% C_dyn_list: list of bearing dynamic load rating [kN]
% C_st_list: list of bearing static load rating [kN]
% f0_list: list of bearing f0 list [-]
% desg_list: list of bearing designations, input is 1x1 cell format from csv

    % values from Fig 11-24 pg 705 machine design:
    % AND from table 9 pg 257 of SKF bearings catalogue
    V = 1.0; % rotation factor, rotating inner ring
    X = 0.56; % radial factor for all deep groove bearings
    Y_list = [2.3,1.99,1.71,1.55,1.45,1.31,1.15,1.04,1.00]; % table 9 pg 257 skf datasheet
    F0Fa_C0_list = [0.172,0.345,0.689,1.03,1.38,2.07,3.45,5.17,6.89];
    e_list = [0.19, 0.22, 0.26, 0.28, 0.30, 0.34, 0.38, 0.42, 0.44];
    regex_seal = "-RS1|-RS2|-RSH|-RSH2|RSJEM"; % bearings with contact seal on one side, pg 258 skf datasheet
    a_skf = 1; % life modification factor, pg 92-98 skf datasheet
    
    bearing_index = -1;
    lifetime = -1;
    for i = 1:length(d_list)
        name_cell = desg_list(i);
        % check if bearing fits case:
        % (dosen"t check for limiting RPM as it"s much higher than current case)
        if d_list(i) >= d_min && ...
           ~isempty(regexp(string(name_cell{1}),regex_seal,"once"))
            e_check_val = (f0_list(i) * abs(F_a)) / (C_0_list(i)*1e3); % tab 9 pg 257 skf datasheet
            [~,Y_index] = closest(F0Fa_C0_list,e_check_val);
            Y = Y_list(Y_index);
            e = e_list(Y_index);
    
            if (abs(F_a)/(F_r*V)) <= e % axal load is irrelevant, pg 256 skf datasheet
                P = F_r; % [N] equivalent load
            else
                P = X * V * F_r + Y * abs(F_a); % [N] equivalent load, eq 11.22a pg 704 machine design
            end
            if P == 0
                warning("Divide by zero")
            end
            % calculate lifetime
            L_P = (a_skf* (K_R*((C_dyn_list(i)*1e3)/P)^3) ) * 1e6; % eq 11.20a pg 701 machine design

            if (L_P > cycles) && (P < C_0_list(i)*1e3)
                % P = P*1e-3
                % (C_dyn_list(i)*1e3)/(P*1e3)
                bearing_index = i; % smallest fitting bearing found, exit loop
                lifetime = L_P;
                break
            end
        end
    end
    if bearing_index == -1
        error("No suitable bearing found")
    end
end