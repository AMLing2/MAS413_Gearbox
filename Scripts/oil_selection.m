clc; close all; clear;

%Parameters to be changed
F_t =   1   ; %Tannhjulskraft             [N]
b =      1  ; %Tannhjulsbredde            [mm]
d1 =     1;  %Delingssirkelens diameter  [mm]
U =     1 ;   %Utveksling                 [-]
v =     1;    %Delingssirkelens hastighet [m/s]


Belastning_hastighetsfaktorer_g1 = ks(F_t, b, d1,U) / v;%[ (N*s) /  (mm^2*m) ]



function ks = ks(F_t, b, d1, U)
    % Beregner k_s til Belastning-hastighetsfaktorer
    % Input:
    %   Ft - Tannhjulskraft (N)
    %   b  - Tannhjulsbredde (mm)
    %   d1 - Delingssirkelens diameter (mm)
    %   U  - Utveksling

    % Output:
    %   ks - (N/mm^2)
    
    ks = (F_t / (b * d1)) * ((U + 1) / U) * 3;
end