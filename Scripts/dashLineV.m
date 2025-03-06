function dashLineV(L, figNr, spRow, spCol)
    % Add dashed line to all figures
    for i = 1 : figNr
        if i < figNr
            figure(i)
            for row = 1 : spRow
                for col = 1 : spCol
                    spIdx = (row-1) * spCol + col;
                    subplot(spRow,spCol,spIdx)
                    xline(L, 'm--', 'LineWidth', 1)
                end
            end
        else
            figure(i)
            xline(L, 'm--', 'LineWidth', 1)
        end
    end
end