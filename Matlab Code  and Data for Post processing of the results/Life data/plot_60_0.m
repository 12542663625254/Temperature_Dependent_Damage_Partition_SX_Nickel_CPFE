clc; clear; close all;

%% 1. EXPERIMENTAL DATA
strain_760 = [0.024, 0.020, 0.016, 0.014]; 
life_760   = [25, 74, 721, 6840];           

strain_980 = [0.024, 0.020, 0.016, 0.014]; 
life_980   = [45, 65, 407, 2225];           

%% 2. SECOND-ORDER POLYNOMIAL FIT (Log-Log Space)
% 760C: Grid starts exactly at first data point life
fit_grid_760 = logspace(log10(life_760(1)), log10(life_760(end)), 200);
p760 = polyfit(log10(life_760), log10(strain_760), 2);
trend_760 = 10.^polyval(p760, log10(fit_grid_760));

% 980C: Start grid at 35 (instead of 45) to ensure curve reaches 0.024 height
fit_grid_980 = logspace(log10(35), log10(life_980(end)), 200); 
p980 = polyfit(log10(life_980), log10(strain_980), 2);
trend_980 = 10.^polyval(p980, log10(fit_grid_980));

%% 3. VISUALIZATION
figure('Color', 'w', 'Position', [200, 200, 800, 700]);
hold on;

% Plot Trendlines
plot(fit_grid_760, trend_760, 'k-', 'LineWidth', 2.5);
plot(fit_grid_980, trend_980, 'r-', 'LineWidth', 2.5);

% Plot Data Points
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
y_labels = arrayfun(@(x) sprintf('%.3f', x), y_ticks, 'UniformOutput', false);
set(gca, 'YTickLabel', y_labels);

% Set limits exactly to match your requirement
xlim([10, 10000]);
ylim([0.014, 0.024]); 

% Labels
xlabel('Cycles to Failure, N_f', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Strain range, \Delta\epsilon', 'FontSize', 26, 'FontName', 'Times New Roman');

%% --- LEGEND ---
h1_leg = plot(nan, nan, 'ks', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
h2_leg = plot(nan, nan, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
lgd = legend([h1_leg, h2_leg], {'760^\circC', '980^\circC'}, ...
             'Location', 'southwest', ...
             'FontSize', 22, ...
             'FontName', 'Times New Roman');
title(lgd, 'DD6, R_{\epsilon}= -1, (60/0)', 'FontName', 'Times New Roman','FontWeight', 'normal');
set(lgd, 'EdgeColor', 'k', 'LineWidth', 1.2);

axis square;
grid off;
hold off;

fprintf('Updated: 980C trendline now extends to touch the 0.024 marker.\n');

%% 9. SAVE FIGURES AS PDFs
% Save Figure 1 (Stress-Strain Hysteresis)
fig1 = figure(1); 
exportgraphics(fig1, 'plot_60_0.pdf', 'ContentType', 'vector');