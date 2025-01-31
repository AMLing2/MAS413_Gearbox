clc;clear;close all
[z_1f, z_2f, z_3f, z_4f, i_totf] = grat2stage(4.5,5,0.001,20,17.3);
z_1f
z_2f
z_3f
z_4f
i_totf

% some outputs: (z = number of teeth)
% z_1 = gear 1 (stage 1), z_2 = gear 2 (stage 1)
% z_3 = gear 3 (stage 2), z_4 = gear 4 (stage 2)

% grat2stage(4.5,5,0.001,18,17.3);
% z_1: 18 , z_2: 79 
% z_3: 18 , z_4: 71 
% i_tot: 17.311728

% grat2stage(4.5,5,0.001,19,17.3); %closest
% z_1: 19 , z_2: 90 
% z_3: 23 , z_4: 84 
% i_tot: 17.299771

% grat2stage(4.5,5,0.001,20,17.3);
% z_1: 20 , z_2: 97 
% z_3: 30 , z_4: 107 
% i_tot: 17.298333 









