close all; clear; clc;


% TO DO:
% - Complete check if stress due to shear can be necglected


% Common input parameters (for all shafts)
n_f = 2; % Safety factor
material = 4140; % (1045 4130 4140 4340)
load_type = "Complex axial";  % ("Pure bending" "Pure axial" "Pure torsion" "Complex axial" "Complex non axial");
surface_finish = "Machined"; % ("Ground" "Machined" "Hot-rolled" "As-forged") For other types, see Machine Design page 368, Figure 6-26
reliability = 95; % [%] reliability factor (50 90 95 99 99.9 99.99 99.999 99.9999)
operating_temperature = 70; % Celsius, defined by Kjell (only significant if > 450)

% Cross sections,  cs_ = [diameter (mm), fillet radius (mm), D_d, P_x (N), M_z (Nmm), M_y (Nmm), T (Nmm)]

%%%%% Shaft 1 %%%%%  From loadingDiagrams_shaft1.m
% Input parameters
r_S1 = 2; % [mm]
D_d_1 = 1.2; % !! ASK ABOUT THIS !!

cs_A =  [30 0.25 D_d_1 0 0 0 83*1e3];
cs_B =  [29 r_S1 D_d_1 0 0 0 cs_A(7)];
cs_G1 = [35 r_S1 D_d_1 -922 66*1e3 -29*1e3 cs_A(7)];
cs_C =  [24 r_S1 D_d_1 cs_G1(2) 0 0 0];

%%%%% Shaft 2 %%%%%  From loadingDiagrams_shaft2.m
% Input parameters
r_S2 = 2; % [mm]
D_d_2 = 1.2; % !! ASK ABOUT THIS !!

cs_E =  [20 r_S2 D_d_2 1757 0 0 0];
cs_G3 = [20 r_S2 D_d_2 cs_E(4) -375*1e3 -128*1e3 362];
cs_G2 = [20 r_S2 D_d_2 -922 -232*1e3 29*1e3 cs_G3(7)];
cs_D =  [20 r_S2 D_d_2 0 0 0];

%%%%% Shaft 3 %%%%%  From loadingDiagrams_shaft3.m
% Input parameters
r_S3 = 2; % [mm]
D_d_3 = 1.2; % !! ASK ABOUT THIS !!

cs_F =  [20 r_S3 D_d_3 -46359 0 0 0];
cs_G4 = [20 r_S3 D_d_3 cs_F(4) 5325*1e3 4169*1e3 24672];
cs_G =  [20 r_S3 D_d_3 0 0 0 cs_G4(7)];
cs_H =  [20 0.4 D_d_3 0 0 0 cs_G4(7)];


%%%%%%%%%%%%%%%% Select corss section to evaluate %%%%%%%%%%%%%%%%
cs = cs_A;
first_itteration = "n"; % ("y" / "n") First itteration for diameter equation (limited geomertry data)
keyseat = "y";          % ("y" / "n") Stress concentration for cross section A and H 

K_t_keyseat = 3.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run("fatigue.m")

%%%%% Modified Goodman Diagram %%%%%
modifiedGoodman(S_y, S_yc, S_ut, S_e, sigma_vm_mean, sigma_vm_amp);