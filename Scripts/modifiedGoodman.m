function modifiedGoodman(S_y, S_yc, S_ut, S_e, sigma_vm_mean, sigma_vm_amp)
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
    
    % Static safety factor
    if sigma_vm_mean >= 0
        n_y = S_y/(sigma_vm_mean + sigma_vm_amp);
    else
        n_y = S_y/abs(sigma_vm_mean - sigma_vm_amp);
    end

    if n_y > 1
        fprintf('n_y = %.2f --> No static failure\n', n_y)
    else
        fprintf('n_y = %.2f --> Static failure\n', n_y)
    end
    
    % Fatigue safety factor
    if sigma_vm_mean >= 0 && sigma_vm_mean + sigma_vm_amp < S_y
        sigma_rev = sigma_vm_amp/(1-(sigma_vm_mean/S_ut));
    elseif sigma_vm_mean >= 0 && sigma_vm_mean + sigma_vm_amp > S_y
        sigma_rev = S_y;
    elseif sigma_vm_mean < 0 && abs(sigma_vm_mean - sigma_vm_amp) < abs(S_yc)
        sigma_rev = sigma_vm_amp;
    elseif sigma_vm_mean < 0 && abs(sigma_vm_mean - sigma_vm_amp) > abs(S_yc)
        sigma_rev = abs(S_yc);
    else
        fprintf('\sigma_{rev} error in Modified Goodman Diagram')
    end
    n_f = S_e/sigma_rev;

    if n_f > 1
        fprintf('n_f = %.2f --> No fatigue failure, infinete life\n', n_f)
    else
        fprintf('n_f = %.2f --> Fatigue failure, finete life\n', n_f)
    end
    
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

    % Calculated point
    plot(sigma_vm_mean, sigma_vm_amp, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
    
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
