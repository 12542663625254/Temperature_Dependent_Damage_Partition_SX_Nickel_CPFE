clc; clear; close all;

%%%%%%% Figure-13 (a)


%% 1. EXPERIMENTAL DATA (30/30 Configuration)
% Replace these arrays with your actual 30/30 experimental data
strain_760 = [0.024, 0.020, 0.016, 0.014]; 
life_760   = [10, 70, 462, 1300];           

strain_980 = [0.024, 0.020, 0.016, 0.014]; 
life_980   = [16, 50, 312, 1407];           

%% 2. SECOND-ORDER POLYNOMIAL FIT (Log-Log Space)
fit_grid = logspace(log10(10), log10(20000), 100);

p760 = polyfit(log10(life_760), log10(strain_760), 2);
trend_760 = 10.^polyval(p760, log10(fit_grid));

p980 = polyfit(log10(life_980), log10(strain_980), 2);
trend_980 = 10.^polyval(p980, log10(fit_grid));

%% 3. VISUALIZATION
figure('Color', 'w', 'Position', [200, 200, 800, 700]);
hold on;

% Plot Trendlines
plot(fit_grid, trend_760, 'k-', 'LineWidth', 2.5);
plot(fit_grid, trend_980, 'r-', 'LineWidth', 2.5);

% Plot Data Points (Large markers for visibility)
scatter(life_760, strain_760, 250, 'ks', 'filled', 'MarkerEdgeColor', 'k');
scatter(life_980, strain_980, 250, 'ro', 'filled', 'MarkerEdgeColor', 'r');

% --- AXIS FORMATTING (Times New Roman) ---
set(gca, 'XScale', 'log', 'YScale', 'log', ...
         'FontSize', 22, ...
         'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, ...
         'TickDir', 'in', ...
         'Box', 'on', ...
         'XMinorTick', 'on', ...
         'YMinorTick', 'on');

% Detailed Y-Axis Ticks
y_ticks = 0.010:0.002:0.028; 
set(gca, 'YTick', y_ticks);

% Formatting labels to 3 decimal places
y_labels = arrayfun(@(x) sprintf('%.3f', x), y_ticks, 'UniformOutput', false);
set(gca, 'YTickLabel', y_labels);

% Limits (Adjusted for better data framing)
xlim([5, 10000]);
ylim([0.014, 0.024]);

% Labels
xlabel('Cycles to Failure, N_f', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Strain range, \Delta\epsilon', 'FontSize', 26, 'FontName', 'Times New Roman');

%% --- LEGEND (30/30 CONFIGURATION) ---
% Dummy handles for legend control
h1_leg = plot(nan, nan, 'ks', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
h2_leg = plot(nan, nan, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r');

lgd = legend([h1_leg, h2_leg], {'760^\circC', '980^\circC'}, ...
             'Location', 'southwest', ...
             'FontSize', 22, ...
             'FontName', 'Times New Roman','FontWeight', 'normal');

% Title updated to reflect 30/30 configuration
title(lgd, 'DD6, R_{\epsilon}= -1, (30/30)', 'FontName', 'Times New Roman','FontWeight', 'normal');
set(lgd, 'EdgeColor', 'k', 'LineWidth', 1.2);

axis square;
grid off;
hold off;

fprintf('30/30 Plot generated with Times New Roman formatting.\n');

%% 9. SAVE FIGURES AS PDFs
% Save Figure 1 (Stress-Strain Hysteresis)
fig1 = figure(1); 
exportgraphics(fig1, 'plot_30_30.pdf', 'ContentType', 'vector');