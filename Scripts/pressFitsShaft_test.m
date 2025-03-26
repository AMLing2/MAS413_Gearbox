clc;close all;clear;
d_hub_outer = 41; % [mm]
d_shaft = 40; % [mm]
l_hub = 37; % [mm]
mu = 0.74; % static dry, mild steel on mild steel, tab 7-1 pg 464 machine design
E_mat = 200; % [GPa]
V_mat = 0.29;
[p,T_max,d_c,h,s,heat_temp_hub,cool_temp_shaft,sigma_t_s,sigma_t_o,sigma_r_s,sigma_r_o] = ...
    pressFitsShaft(d_hub_outer,d_shaft,l_hub,mu,E_mat,E_mat,V_mat,V_mat);

bearing_inner = 70; % [mm]
[p1,T_max1,d_c1,h1,s1,heat_temp1,cool_temp1,sigma_t_s1,sigma_t_o1,sigma_r_s1,sigma_r_o1] = ...
    pressFitsHub(d_hub_outer,d_shaft,l_hub,mu,E_mat,E_mat,V_mat,V_mat,bearing_inner);

datatable = table([p,T_max,d_c,h,s,heat_temp_hub,cool_temp_shaft,sigma_t_s,sigma_t_o,sigma_r_s,sigma_r_o]' ...
    , ["p","T_max","d_c","h","s","heat_temp","cool_temp_shaft","sigma_t_s","sigma_t_o","sigma_r_s","sigma_r_o"]', 'VariableNames', {'data', 'name'});
disp(datatable);

datatable1 = table([p1,T_max1,d_c1,h1,s1,heat_temp1,cool_temp1,sigma_t_s1,sigma_t_o1,sigma_r_s1,sigma_r_o1]' ...
    , ["p1","T_max1","d_c1","h1","s1","heat_temp1","cool_temp1","sigma_t_s1","sigma_t_o1","sigma_r_s1","sigma_r_o1"]', 'VariableNames', {'data', 'name'});
disp(datatable1);