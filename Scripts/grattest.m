clc;clear;close all
[z_1f, z_2f, z_3f, z_4f, i_totf] = grat2stage(4,5,0.001,19,13.7);
z_1f
z_2f
z_3f
z_4f
i_totf

% some outputs: (z = number of teeth)
% z_1 = gear 1 (stage 1), z_2 = gear 2 (stage 1)
% z_3 = gear 3 (stage 2), z_4 = gear 4 (stage 2)

% grat2stage(4,5,0.001,18,13.7)
% z_1f = 18
% z_2f = 73
% z_3f = 29
% z_4f = 98
% i_totf = 13.704981

% grat2stage(4,5,0.001,19,13.7) % second closest
% z_1f = 19
% z_2f = 83
% z_3f = 22
% z_4f = 69
% i_totf = 13.700957

% grat2stage(4,5,0.001,20,13.7)
% z_1f = 20
% z_2f = 87
% z_3f = 20
% z_4f = 63
% i_totf = 13.702500

%grat2stage(4,5,0.001,21,13.7) % closest to 13.7 but z_1 > 20
% z_1f = 21
% z_2f = 88
% z_3f = 26
% z_4f = 85
% i_totf = 13.699634















% z_1 = 18
% z_2 = 79
% z_3 = 25
% z_4 = 78
% i_tot = 13.693333

% z_1 = 19
% z_2 = 94
% z_3 = 13 %too low....
% z_4 = 36
% i_tot: 13.700405

% z_1 = 20
% z_2 = 87
% z_3 = 20
% z_4 = 63
% i_tot = 13.702500

% z_1 = 21 % closest one to 13.7 but z_1 > 20
% z_2 = 88
% z_3 = 26
% z_4 = 85
%i_totf = 13.699634
