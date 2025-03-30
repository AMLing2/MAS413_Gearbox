%%% Outputs
% L_sh : length of shoulder [mm]
%%% Inputs
% a : axial force on the shoulder [N]
% d_shaft : diameter of the shaft before the shoulder [mm]
% d_shoulder : diameter of the shoulder [mm]
% n_f : saftey factor [mm]
% S_y : Yield strength of material [MPa]

function L_sh = shoulder_length(a,d_shaft,d_shoulder,n_f,S_y)
    A_norm = pi * (d_shoulder/2)^2 - pi * (d_shaft/2)^2;
    sigma_x = a / (A_norm);
    L_sh = a/ (sqrt( (((S_y/n_f) *sqrt(2))^2 - 2 * sigma_x^2) / (6)) * 2 * pi * ((d_shoulder - d_shaft)/2) );
end