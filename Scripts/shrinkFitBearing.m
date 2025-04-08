% chapter 10.12 machine design pg 620
% case: known hub inner diameter (bearing selected), increase shaft outer diamter to fit the hub
% (for bearings)
%%%INPUTS:
% d_h_o : outer diameter of hub / gear (D for bearings) [mm]
% d_h_i : inner diameter of hub / bearing (d for bearings) [mm]

%%%OUTPUTS:
% d_s_c : changed diameter of shaft [mm]
% h : gear hole diameter tolerence
% s : shaft diameter tolerence
% heat_temp_hub : temperature of hub / bearing needed for the fit [deg C]
% temp_shaft : temperature of shaft needed for the fit [deg C]
% bearing_clearance : radial clearance in the bearing (likely incorrect) [um]

function [d_s_c,h,s,heat_temp_bearing,temp_shaft,bearing_clearance] ...
    = shrinkFitBearing(d_h_o,d_h_i,fit)

    % class 8 fit (h8), Appendix E-1 fundamentals of machine component design pg 854
    C_h = 0.0052;
    C_s = 0.0052;
    C_i_8 = 0.0010;
    C_i_7 = 0.0005;
    % find values for a k5 fit
    if strcmp(fit,"h7s6")
        i = C_i_8 * d_h_i; % [mm] average interference class 8, h7s6
        % h7s6 from table 4.1a (Steinschaden Lecture 6 slide 12 of
        % "078-Engineering-Drawings-Lecture-Linear-Fits-Tolerances.pdf")
    elseif strcmp(fit,"h7p6")
        i = C_i_7 * d_h_i; % [mm] average interference class 7, h7p6
    else
        error("select fit h7s6 or h7p6")
    end
    % i_7 = C_i_7 * d_h_i; % [mm] average interference h7
    % i_8 = C_i_8 * d_h_i; % [mm] average interference h8
    
    d_s_c = d_h_i + i; % [mm] new radius of shaft
    h = C_h * d_h_i^(1/3); % gear hole diameter tolerence 
    s = C_s * d_s_c^(1/3); % shaft diameter tolerence 

    % temperature calculations: pg 577 university physics book
    temp_room = 22; % room temperature [deg C]
    temp_shaft = temp_room;
    L_0 = 2*pi*(d_h_i/2); % initial circumference [mm]
    fit_tol = s+h; % extra radius for easier fit [mm]
    L_heated = 2*pi*((d_h_i + i + fit_tol)/2); % finishing length for fit [mm]
    alpha = 1.2e-5; % Coefficient of liear expansion [1/deg C], tab 17.1 pg 578
    heat_temp_bearing = ((L_heated - L_0) / (alpha * L_0)) + temp_room; % [deg C] eq 17.6 pg 576

    max_heating = 110; % bad for bearings to be above 125 deg C
    if heat_temp_bearing > max_heating
        %warning("bearing d = %dmm has to be cooled for interference fit",d_s)
        temp_shaft = temp_room - (heat_temp_bearing - max_heating); % cool down shaft instead
        heat_temp_bearing = max_heating;
    end

    % bearing clearance reduction calcs
    temp_outer = 68; % [deg C] temperature of outer ring, estimate
    working_temp = 70; % [deg C]
    d_m = (d_h_o + d_h_i) / 2; % [mm] bearing mean diameter
    delta_r_temp = 0.012 * (working_temp - temp_outer) * d_m; % [um], converts to [um] in the datasheet
    f1 = 0.4773 * (d_h_i/d_h_o) + 0.5507; % reduction factor for inner ring, diagram 2 pg184 skf catalog assuming linear func
    delta_r_fit = f1 * i * 1e3; % [um] assuming loose fit on outer ring seat
    r_op = 0; % [um] required operating clearance, equal to 0 for ball bearings
    r = r_op + delta_r_temp + delta_r_fit; % pg 184 skf datasheet
    
    clearance_list = [2 2 3 5 5 6 8 10 12 15 18 18 20 25]; % pg 252 table 6 skf catalog
    bore_diam_list = [2.5 6 10 18 24 30 40 50 65 80 100 120 140 160 180 200];
    [~,i] = closest(bore_diam_list,d_h_i);
    min_clearance = clearance_list(i);
    bearing_clearance = min_clearance - r; % [um]
    % if bearing_clearance < 0
    %      warning("bad bearing clearance")
    % end
end