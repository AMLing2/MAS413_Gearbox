function plotLD(x, f, colFill)
    % Plot f(x) as a black line with shaded area 
    % between f(x) and abscissa

    % Aesthetics :sparkly_glitter:
    colLine = 'k';
    lineW = 2;

    % Filled Area
    X = [x, fliplr(x)];
    Y = [f, zeros(size(f))];

    % Actual Plotting
    hold on
    fill(X, Y, colFill, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(x, f,'Color', colLine','LineWidth', lineW)
    yline(0, 'k')

    % Detect Jumps and Draw Line
    jumpTol = 10;
    jumps = find( abs( diff(f) ) > jumpTol );
    for i = 1 : length(jumps)
        xJump = x(jumps(i) + 1);
        plot( [ xJump xJump ], [ f(jumps(i), f(jumps(i)+1) )] )
    end

end