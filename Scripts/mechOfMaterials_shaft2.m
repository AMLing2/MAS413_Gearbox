 clc; clear; close all;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % MAS413 Project: Mechanics of Materials - Shaft 2%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Common Plotting Constants
 colFill = [0.7765 0.9176 0.9843];
 resolution = 100;
 wPlot = 22;
 hPlot = 16;

 % Given information from project description - Needs adjustment of values
 i_tot = 17.3; % [-]
 i_1 = 1; % [-]
 i_2 = 2; % [-]


 n_in = 1450; % [RPM]
 omega_in = n_in * 2*pi / 60; % [rad/sec]
 n_out = (n_in/i_tot); %[RPM]
 omega_out = n_out * 2*pi / 60; %[rad/sec]
 
 eta = 0.96; %[-] Efficiency
 eta_tot = eta^2;

 P_in = 12.5e3; % [W]
 P_out = P_in*eta_tot; %[W] (Squared efficiency because there are two stages)
 T_M = P_in/omega_in; %[Nm]
 T_out = P_out/omega_out; %[Nm] 
 
 alpha = 20; % [degrees] Helix Angle
 beta = 15;  % [degrees] Pressure Angle

 % Lengths on shaft - Needs adjustment of value
 L_EG3  = 0.05; % [m]
 L_G3G2 = 0.10; % [m]
 L_G2D = 0.15; % [m]
 L_tot = L_EG3 + L_G3G2 + L_G2D; % [m]
 L_G3D = L_G3G2+L_G2D; % [m]
 
 % Radius of gears - Calculated Elsewhere - Needs adjustment of value
 r_G1 = 0.25; % [m]
 r_G2 = 0.25; % [m]
 r_G3 = 0.25; % [m]
 r_G4 = 0.25; % [m]

 %% Calculations
 % Calculated lengths & Torques
 L_E_G2 = L_EG3 + L_G3G2; % [m]
 
 %Gear 2 forces
 F_t2 = T_M/r_G1; % [N]
 F_a2 = F_t2 * tand(beta); % [N]
 F_r2 = F_t2 * tand(alpha)/cosd(beta); % [N]
 
 % Gear 3 forces
 F_t3 = (T_out*i_tot) / r_G4; % [N]
 F_a3 = F_t3 * tand(beta); % [N]
 F_r3 = F_t3 * tand(alpha)/cosd(beta); % [N]
 
 %% Reaction forces @ bearings
 L_ED = L_tot;

 F_Ex = F_a2-F_a3;

 F_Ez = (F_r3*L_G3D - F_r2*L_G2D - F_a3*r_G3 - F_a2*r_G2)/L_ED;

 F_Ey= (F_t3*L_G3D + F_t2*L_G2D)/L_ED;

 
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
 
  x = linspace(0, L_EG3, resolution);
  P = ones(size(x)) * F_Ex; % [N]
  V = ones(size(x)) * (-F_Ey); % [N]
  M = ones(size(x)) .* (-F_Ey * x); % [Nm]
  T = zeros(size(x)); % [Nm]
 
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
 
 x = linspace(L_EG3, L_G3G2, resolution);
 P = ones(size(x))* -(F_Ex+F_a3); % [N]
 V = ones(size(x)) * -(F_Ey+F_t3); % [N]
 M = -F_Ey*x + F_t3*(x - L_EG3); % [Nm]
 T = ones(size(x)) * F_t3*r_G3; % [Nm]
 
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
 
 x = linspace(L_G3G2, L_G2D, resolution);
 P = ones(size(x)) * (F_Ex - F_a2 - F_a3); % [N]
 V = ones(size(x)) * (F_t2 - F_Ey + F_t3); % [N]
 M = F_t2 * (x- L_E_G2) - F_Ey*x + F_t3 * (x - L_EG3); % [Nm]
 %T = ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2); % [Nm]
 
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
 
% %% XZ - Plane
 
 % Figure setup
 resolution = 100;
 figHandle = 2;
 xPos = 10;
 yPos = 3;
 wPlot = 22;
 hPlot = 16;
 
 XZplaneFig = figure(figHandle);
 set(figHandle,'Units','Centimeter')
 set(figHandle,'Position',[xPos yPos wPlot hPlot]);
 sgtitle('XZ - Plane')
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
 
 x = linspace(0, L_EG3, resolution);
 P = ones(size(x))*(-F_Ex); % [N]
 V = ones(size(x))*(-F_Ez); % [N]
 M = ones(size(x)).* (x *(-F_Ex)); % [Nm]
 T = zeros(size(x)); % [Nm]
 
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
 x = linspace(L_EG3, L_G3G2, resolution);
 P = ones(size(x))* -(F_Ex+F_a3); % [N]
 V = ones(size(x)) * (F_r3-F_Ez); % [N]
 M = F_Ez*x - F_r3*(x - L_EG3); % [Nm]
 T = ones(size(x)) * F_t3*r_G3; % [Nm]
 
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
 
 x = linspace(L_G3G2, L_G2D, resolution);
 P = ones(size(x)) * (F_a2 - F_a3 - F_Ex); % [N]
 V = ones(size(x)) * (F_r3 - F_Ez - F_r2); % [N]
 M = F_Ez * x - F_r3*(x - L_EG3) + F_a3 * r_G3 * (x - L_EG3) + F_a2*r_G2*(x-L_G3G2); % [Nm]
 T = ones(size(x)) * (F_t3*r_G3 - F_t2*r_G2); % [Nm]
 
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