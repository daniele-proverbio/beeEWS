% Plot PSD and fit 
function plot_psd_fit(nc, label)
    log_f   = log10(nc.freq);
    log_p   = log10(nc.psd);
    fit_line = nc.fit_slope * log_f + nc.fit_offset;

    figure('Name','PSD Noise Color Fit','Color','w','Position',[100 100 700 420]);

    % Raw PSD
    semilogy(nc.freq, nc.psd, 'Color',[0.6 0.6 0.6], 'LineWidth',1);
    hold on;

    % Fitted line (back in linear scale)
    semilogy(nc.freq, 10.^fit_line, 'r-', 'LineWidth', 2.5);

    % Highlight fit region
    semilogy(nc.freq(nc.idx_fit), nc.psd(nc.idx_fit), ...
             'b-', 'LineWidth', 1.4);

    set(gca,'XScale','log');
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density');
    title(sprintf('%s  |  Noise: %s  |  \\beta = %.2f  |  R² = %.3f', ...
                  label,nc.color, nc.beta, nc.r_squared), 'Interpreter','tex');
    legend('Full PSD','Power-law fit (1/f^{\beta})','Fit region','Location','southwest');
    hold off;
end