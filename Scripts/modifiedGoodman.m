function modifiedGoodman(S_yc, S_ut, S_e)
    % Modified-Goodman Graph (Maskinelementer, lecture 4 slide 12-22)
    % (Work in progress)

    % Define mean stress range
    sigma_mean_goodman_L = linspace(S_yc, 0, 100); % Left side (compressive)
    sigma_mean_goodman_R = linspace(0, S_ut, 100); % Right side (tensile)
    
    % Modified-Goodman equation
    sigma_amp_goodman_L = S_e * ones(size(sigma_mean_goodman_L));
    sigma_amp_goodman_R = S_e * (1 - sigma_mean_goodman_R / S_ut);
    
    % Static yeilding line (first cycle)
    S_y_goodman_L = S_y * (1 - sigma_mean_goodman_L / S_yc);
    S_y_goodman_R = S_y * (1 - sigma_mean_goodman_R / S_y);
    
    % Check for static failure
    % n_y = S_y/simga_max;
    
    % Check for fatigue failure
    % sigma_rev = sigma_amp_vm/(1-(sigma_mean_vm/S_ut));
    % n_f = S_e/simga_rev;
    
    % Intersecting points for indexing
    [~, intersect_L] = min(abs(sigma_amp_goodman_L - S_y_goodman_L));
    [~, intersect_R] = min(abs(sigma_amp_goodman_R - S_y_goodman_R));
    
    figure; hold on;
    % Colour area between graphs and x-axis
    fill([sigma_mean_goodman_L, sigma_mean_goodman_R], [S_y_goodman_L(1:intersect_L-1), sigma_amp_goodman_L(intersect_L:end),...
        sigma_amp_goodman_R(1:intersect_R-1), S_y_goodman_R(intersect_R:end)], [0.3010 0.7450 0.9330], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
    
    % Equations
    plot(sigma_mean_goodman_L, sigma_amp_goodman_L, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2); hold on
    plot(sigma_mean_goodman_R, sigma_amp_goodman_R, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
    plot(sigma_mean_goodman_L, S_y_goodman_L, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
    plot(sigma_mean_goodman_R, S_y_goodman_R, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
    
    % Points for strengths
    plot(0, S_e, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_e point
    plot(S_ut, 0, 'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerSize', 6); % S_ut point
    plot(0, S_y, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
    plot(S_y, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_y point
    plot(S_yc, 0, 'ko', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerSize', 6); % S_yc point
    axis([S_yc-100 S_ut+100 0 S_y+100]);
    
    % Axis labels ! NEEDS ADJUSTMENT
    xticks([S_yc, S_y, S_ut]); % Set x-axis tick positions
    yticks([S_e, S_y]); % Set y-axis tick positions
    xticklabels({'S_{yc}', 'S_{y}', 'S_{ut}'}); % Custom x-axis labels
    yticklabels({'S_e', 'S_y'}); % Custom y-axis labels
    ax = gca;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    ax.YAxisLocation = 'origin'; % Move y-axis to x = 0
    
    xlabel('\sigma_m [MPa]', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('\sigma_a [MPa]', 'FontSize', 14, 'FontWeight', 'bold');
end
