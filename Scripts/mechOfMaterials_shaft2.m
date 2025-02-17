 clc; clear; close all;
 % MAS413 Project: Mechanics of Materials - Shaft 2
 
 %% Constants
 
 % Given information from project description - Needs adjustment of value
 n_1 = 1450; % [RPM]
 omega_1 = n_1 * 2*pi / 60; % [rad/sec]
 P_1 = 12.5e3; % [W]
 i_tot = 17.3; % [-]
 i_1 = 1; % [-]
 i_2 = 2; % [-]
 alpha = 20; % [degrees] Helix Angle
 beta = 15;  % [degrees] Pressure Angle
 


 % Lengths on shaft - Needs adjustment of value
 L_E_G3  = 0.05; % [m]
 L_G3_G2 = 0.10; % [m]
 L_G2_D = 0.15; % [m]
 L_tot = L_E_G3 + L_G3_G2 + L_G2_D; % [m]
 
 % Radius of gears - Calculated Elsewhere - Needs adjustment of value
 r_G2 = 0.25  % [m]
 r_G3 = 0.25; % [m]
 
 
 % Calculated lengths & Torques
 L_E_G2 = L_E_G3 + L_G3_G2; % [m]
 
 T_M = P_1 / omega_1; % [Nm]
 T_S2 = T_M * i1; %[Nm]
 
 %Gear 2 forces
 F_t2 = T_S2/r_G2; % [N]
 F_a2 = F_t2 * tand(beta); % [N]
 F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]
 
 % Gear 3 forces
 F_t3 = T_S2 / r_G3; % [N]
 F_a3 = F_t3 * tand(beta); % [N]
 F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]
 
 % Reaction forces @ bearings
 L_BC = L_BG1 + L_G1C; % [m]

 F_By = F_t1*L_G1C/L_BC; % [N]
 F_Bz = F_r1*L_G1C/L_BC; % [N]
 
 
 %% XY - Plane
 
 % Figure setup
 resolution = 100;
 figHandle = 1;
 xPos = 10;
 yPos = 3;
 wPlot = 22;
 hPlot = 16;
 
 XYplaneFig = figure(figHandle);
 set(figHandle,'Units','Centimeter')
 set(figHandle,'Position',[xPos yPos wPlot hPlot]);
 sgtitle('XY - Plane')
 subplot(2,2,1)
 xlabel('[m]')
 ylabel('[N]')
 title('Axial Force P(x)')
 subplot(2,2,2)
 xlabel('[m]')
 ylabel('[N]')
 title('Shear Force V(x)')
 subplot(2,2,3)
 xlabel('[m]')
 ylabel('[Nm]')
 title('Bending Moment M(x)')
 subplot(2,2,4)
 xlabel('[m]')
 ylabel('[Nm]')
 title('Axial Torque T(x)')
 
 % 0 < x < L_E_G3
 
 x = linspace(0, L_AB, resolution);
 %P = zeros(size(x)); % [N]
 %V = zeros(size(x)); % [N]
 %M = zeros(size(x)); % [Nm]
 %T = ones(size(x)) * T_S2; % [Nm]
 
 subplot(2,2,1)
 hold on
 plot(x,P,'r','LineWidth',2)
 subplot(2,2,2)
 hold on
 plot(x,V,'r','LineWidth',2)
 subplot(2,2,3)
 hold on
 plot(x,M,'Color','#FF8800','LineWidth',2)
 subplot(2,2,4)
 hold on
 plot(x,T,'Color','#FF8800','LineWidth',2)
 
 % L_E_G3 < x < L_G3_G2
 
 x = linspace(L_E_G3, L_G3_G2, resolution);
 %P = zeros(size(x)); % [N]
 %V = ones(size(x)) * (-F_By); % [N]
 %M = F_By * (x - L_AB); % [Nm]
 %T = ones(size(x)) * T_S2; % [Nm]
 
 subplot(2,2,1)
 hold on
 plot(x,P,'r','LineWidth',2)
 subplot(2,2,2)
 hold on
 plot(x,V,'r','LineWidth',2)
 subplot(2,2,3)
 hold on
 plot(x,M,'Color','#FF8800','LineWidth',2)
 subplot(2,2,4)
 hold on
 plot(x,T,'Color','#FF8800','LineWidth',2)
 
 % L_G3_G2 < x < L_G2_D
 
% x = linspace(L_G3_G2, L_AC, resolution);
% P = ones(size(x)) * (-F_a1); % [N]
% V = ones(size(x)) * (F_By - F_t1); % [N]
% M = F_By * (x - L_AB) - F_t1 * (x - L_AG1); % [Nm]
% T = ones(size(x)) * (T_S2 - F_t1*r_G3); % [Nm]
% 
% subplot(2,2,1)
% hold on
% plot(x,P,'r','LineWidth',2)
% subplot(2,2,2)
% hold on
% plot(x,V,'r','LineWidth',2)
% subplot(2,2,3)
% hold on
% plot(x,M,'Color','#FF8800','LineWidth',2)
% subplot(2,2,4)
% hold on
% plot(x,T,'Color','#FF8800','LineWidth',2)
% 
% %% XZ - Plane
% 
% % Figure setup
% resolution = 100;
% figHandle = 2;
% xPos = 10;
% yPos = 3;
% wPlot = 22;
% hPlot = 16;
% 
% XZplaneFig = figure(figHandle);
% set(figHandle,'Units','Centimeter')
% set(figHandle,'Position',[xPos yPos wPlot hPlot]);
% sgtitle('XZ - Plane')
% subplot(2,2,1)
% xlabel('[m]')
% ylabel('[N]')
% title('Axial Force P(x)')
% subplot(2,2,2)
% xlabel('[m]')
% ylabel('[N]')
% title('Shear Force V(x)')
% subplot(2,2,3)
% xlabel('[m]')
% ylabel('[Nm]')
% title('Bending Moment M(x)')
% subplot(2,2,4)
% xlabel('[m]')
% ylabel('[Nm]')
% title('Axial Torque T(x)')
% 
% % 0 < x < L_E_G3
% 
% x = linspace(0, L_E_G3, resolution);
% P = zeros(size(x)); % [N]
% V = zeros(size(x)); % [N]
% M = zeros(size(x)); % [Nm]
% T = ones(size(x)) * T_S2; % [Nm]
% 
% subplot(2,2,1)
% hold on
% plot(x,P,'r','LineWidth',2)
% subplot(2,2,2)
% hold on
% plot(x,V,'r','LineWidth',2)
% subplot(2,2,3)
% hold on
% plot(x,M,'Color','#FF8800','LineWidth',2)
% subplot(2,2,4)
% hold on
% plot(x,T,'Color','#FF8800','LineWidth',2)
% 
% % L_E_G3 < x < L_G3_G2
% 
% x = linspace(L_E_G3, L_G3_G2, resolution);
% P = zeros(size(x)); % [N]
% V = ones(size(x)) * (F_Bz); % [N]
% M = - F_Bz * (x - L_AB); % [Nm]
% T = ones(size(x)) * T_S2; % [Nm]
% 
% subplot(2,2,1)
% hold on
% plot(x,P,'r','LineWidth',2)
% subplot(2,2,2)
% hold on
% plot(x,V,'r','LineWidth',2)
% subplot(2,2,3)
% hold on
% plot(x,M,'Color','#FF8800','LineWidth',2)
% subplot(2,2,4)
% hold on
% plot(x,T,'Color','#FF8800','LineWidth',2)
% 
% % L_G3_G2 < x < L_G2_D
% 
% x = linspace(L_G3_G2, L_G2_D, resolution);
% P = ones(size(x)) * (-F_a1); % [N]
% V = ones(size(x)) * (F_Bz - F_r1); % [N]
% M = F_r1 * (x - L_AG1) - F_Bz * (x - L_AB) - F_a1 * r_G3; % [Nm]
% T = ones(size(x)) * (T_S2 - F_t1*r_G3); % [Nm]
% 
% subplot(2,2,1)
% hold on
% plot(x,P,'r','LineWidth',2)
% subplot(2,2,2)
% hold on
% plot(x,V,'r','LineWidth',2)
% subplot(2,2,3)
% hold on
% plot(x,M,'Color','#FF8800','LineWidth',2)
% subplot(2,2,4)
% hold on
% plot(x,T,'Color','#FF8800','LineWidth',2)