clc; close all; clear;
%% Parameters

d_shaft_1 = 1; %[m]
d_shaft_3 = 1; %[m]

T_shaft1 = 1; %[Nm]
T_shaft3 = 1; %[Nm]

S_yield = 100; %[Pa]

%% Calculated values

r_shaft_1 = d_shaft_1/2; %[m]
r_shaft_3 = d_shaft_3/2; %[m]

w_1 = d_shaft_1/4;  %[Nm]
h_1 = d_shaft_1/6;  %[Nm]
l_1 = 1.5*d_shaft_1;%[Nm]

w_3 = d_shaft_3/4;  %[Nm]
h_3 = d_shaft_3/6;  %[Nm]
l_3 = 1.5*d_shaft_3;%[Nm]


%% Check for shear failure

tau_key_1 = (T_shaft1/r_shaft_1)/(w_1*l_1); %[Pa]
n_shear_1 = (0.577*S_yield)/tau_key_1; %[-]

tau_key_3 = (T_shaft3/r_shaft_3)/(w_3*l_3); %[Pa]
n_shear_3 = (0.577*S_yield)/tau_key_3; %[-]


%% Check for compression failure

sigma_key_1 = (T_shaft1/r_shaft_1)/((h_1/2)*l_1); %[Pa]
n_compression_1 = (0.577*S_yield)/sigma_key_1; %[-]

sigma_key_3 = (T_shaft3/r_shaft_3)/((h_3/2)*l_3); %[Pa]
n_compression_3 = (0.577*S_yield)/sigma_key_3; %[-]

