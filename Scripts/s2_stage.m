clc;clear;close all;
n_1 = 1450; % Speed of input shaft 1 [rpm]
P_1 = 12.5e3; % Effect on shaft [W]
i_tot = 17.3; % total gear ratio
alpha = 20; % pressure angle [deg]
beta = 15; % helix angle [deg]

n_out = n_1/i_tot;


% calculations for internal gears
 i_1 = 5; % first gear ratio of 2 stage from table 15-38, Lec2 pg12
 i_2 = i_tot/i_1

z_1_min = 19; % minimum 18-20 from teacher [# teeth]
z_2 = round(i_1 * z_1_min);
% z_2 = 97;
z_1_divisors = alldivisors(z_1_min);
z_2_divisors = alldivisors(z_2);
stage1_divisors = intersect(z_1_divisors,z_2_divisors);

if ~isempty(stage1_divisors)
    disp(["stage 1 divisors: ",stage1_divisors]);
    warning("stage 1 gears are not relative prime") 
else
    disp(["stage 1 divisors: ",stage1_divisors]);
    % disp("\nStage 1 gears ARE relative")
end

i_2 = i_tot/(z_2/z_1_min)

z_4 = 173; % from converting i_2 to rational
z_3 = 50;

z_3 = 20;
z_4 = round(z_3 * i_2) + 1;

z_4_divisors = alldivisors(z_4);
z_3_divisors = alldivisors(z_3);
stage2_divisors = intersect(z_4_divisors,z_3_divisors);

if ~isempty(stage2_divisors)
    disp(["stage 2 divisors:",stage2_divisors]);
    warning("stage 2 gears are not relative prime")
else
    disp(["stage 2 divisors:",stage2_divisors]);
    % disp("\nStage 1 gears ARE relative")
end

final_i_tot = (z_2/z_1_min) * (z_4/z_3)

function d = alldivisors(n)
    % taken from: https://www.mathworks.com/matlabcentral/answers/21542-find-divisors-for-a-given-number#answer_28371
    K = 1:n;
    d = K(rem(n,K)==0);
    d(1) = []; % remove first element of 1
end
