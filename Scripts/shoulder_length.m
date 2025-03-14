function [L_sh] = shoulder_length(a_left,a_right,d_shaft,d_shoulder,n_f,U_t)
    P = a_left - a_right;
    A_norm = pi * (d_shoulder/2)^2 - pi * (d_shaft/2)^2;
    sigma_x = P / (A_norm);
    L_sh = P/ (sqrt( (((U_t/n_f) *sqrt(2))^2 - 2 * sigma_x^2) / (6)) * 2 * pi * ((d_shoulder - d_shaft)/2) );
end