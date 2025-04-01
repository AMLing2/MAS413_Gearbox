% function modifiedGoodman(S_y, S_yc, S_ut, S_e, sigma_vm_mean, sigma_vm_amp)
function modifiedGoodman(S_y, S_ut, shaft_design_results, titleName)
    
    % Modified-Goodman Graph % MAS236 L4 s12-22)
    
    S_e = min(shaft_design_results(:, 1));
    % S_e = shaft_design_results(:, 1);
    sigma_e_m = shaft_design_results(:, 2);
    sigma_e_a = shaft_design_results(:, 3);
    n_y = shaft_design_results(:, 4);
    n_f = shaft_design_results(:, 5);

    % Define mean stress range
    % sigma_mean_goodman_L = linspace(S_yc, 0, 100); % Left side (compressive)
    sigma_mean_goodman_R = linspace(0, S_ut, 100); % Right side (tensile)
    
    % Static yeilding line (first cycle)
    % S_y_goodman_L = S_y * (1 - sigma_mean_goodman_L / S_yc);
    S_y_goodman_R = S_y * (1 - sigma_mean_goodman_R / S_y);

    % Modified-Goodman equation
    sigma_amp_goodman_R = S_e * (1 - sigma_mean_goodman_R / S_ut);
    
    % % Intersecting points for indexing
    % [~, intersect_L] = min(abs(sigma_amp_goodman_L - S_y_goodman_L));
    [~, intersect_R] = min(abs(sigma_amp_goodman_R - S_y_goodman_R));
    
    figure; hold on;

    % Colour area between graphs and x-axis
    % fill([sigma_mean_goodman_L, sigma_mean_goodman_R], [S_y_goodman_L(1:intersect_L-1), sigma_amp_goodman_L(intersect_L:end),...
    %     sigma_amp_goodman_R(1:intersect_R-1), S_y_goodman_R(intersect_R:end)], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    fill([sigma_mean_goodman_R, 0], [sigma_amp_goodman_R(1:intersect_R-1), S_y_goodman_R(intersect_R:end), 0], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;

    % Equations
    
    % Goodman equation
    plot(sigma_mean_goodman_R, sigma_amp_goodman_R, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)

    % plot(sigma_mean_goodman_L, S_y_goodman_L, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
    plot(sigma_mean_goodman_R, S_y_goodman_R, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
    

    % Points for strengths
    plot(S_ut, 0, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_ut point
    plot(0, S_y, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
    plot(S_y, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
    plot(0, S_e, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_e min point
    % plot(S_yc, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_yc point
    
    % List of colors for loop plotting
    marker_colour = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};

    for i = 1:height(shaft_design_results)

    % Modified-Goodman equation
    % sigma_amp_goodman_L = S_e(i) * ones(size(sigma_mean_goodman_L));
    % sigma_amp_goodman_R = S_e(i) * (1 - sigma_mean_goodman_R / S_ut);
    
    % Goodman equation
    % plot(sigma_mean_goodman_L, sigma_amp_goodman_L, '--', 'Color', marker_colour{i}, 'LineWidth', 2)
    % plot(sigma_mean_goodman_R, sigma_amp_goodman_R, '--', 'Color', marker_colour{i}, 'LineWidth', 2)
    
    % S_e points
    % plot(0, S_e(i), 'x', 'Color', marker_colour{i}, 'MarkerSize', 6); % S_e point

    % Calculated point of interest 
    poi(i) = plot(sigma_e_m(i), sigma_e_a(i), 'ko', 'MarkerFaceColor', marker_colour{i}, 'MarkerSize', 6);
    
    end

    % Axis labels ! NEEDS ADJUSTMENT
    axis([0 S_ut+10 0 S_y+10]);
    % axis([S_yc-100 S_ut+100 0 S_y+100]);
    xticks([S_y, S_ut]); % Set x-axis tick positions
    % xticks([S_yc, S_y, S_ut]); % Set x-axis tick positions
    yticks([0, S_e, S_y]); % Set y-axis tick positions
    % yticks([S_e, S_y]);
    xticklabels({sprintf('S_y = %.0f', S_y), sprintf('S_{ut} = %.0f', S_ut)}); % Custom x-axis labels
    yticklabels({0, sprintf('S_e = %.0f', S_e), sprintf('S_y = %.0f', S_y)}); % Custom y-axis labels
    % yticklabels({'S_e', 'S_y'});
    ax = gca;
    set(gca, 'FontSize', 14);
    ax.YAxisLocation = 'origin'; % Move y-axis to x = 0
    
    % Add manual labels at the ends
    xlabel('\sigma_m [MPa]');
    ylabel('\sigma_a [MPa]');
    title(titleName)
    leg = legend('show');

    if titleName == 'Shaft 1'
        legend(poi, {'A', '0', '1', '2'}, 'NumColumns', 2, 'Location', 'northeast')
        title(leg,'Cross Sections')
    elseif titleName == 'Shaft 2'
        legend(poi, {'3', '4', '5', '6'}, 'NumColumns', 2, 'Location', 'northeast')
        title(leg,'Cross Sections')
    else
        legend(poi, {'7', '8', '9', 'H'}, 'NumColumns', 2, 'Location', 'northeast')
        title(leg,'Cross Sections')
    end
end
