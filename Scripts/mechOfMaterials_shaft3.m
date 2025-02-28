clc; clear; close all;
% MAS413 Project: Mechanics of Materials - Shaft 3

%% Constants

% Common Plotting Constants
colFill = [0.7765 0.9176 0.9843];
resolution = 100;
wPlot = 22;
hPlot = 16;

% Given information
n_1 = 1450; % [RPM]
omega_1 = n_1 * 2*pi / 60; % [rad/sec]
P_1 = 12.5e3; % [W]
i_tot = 17.3; % [-]
i_1 = 1; % [-]
i_2 = 2; % [-]
alpha = 20; % [degrees] Helix Angle
beta = 15;  % [degrees] Pressure Angle

