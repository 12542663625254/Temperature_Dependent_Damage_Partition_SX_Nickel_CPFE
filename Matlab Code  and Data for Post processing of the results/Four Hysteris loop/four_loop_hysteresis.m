clc; clear; close all;


%%%%%%%% Figure- 3 %%%%%%%%


%% 1. LOAD DATA
data_60_0 = readmatrix('new_60_Tension.csv'); 
data_30_30 = readmatrix('new_data_simulated.csv');
data_0_60 = readmatrix('Fatigue_0_60.csv');

strain_60_0 = data_60_0(:, 25);   stress_60_0 = data_60_0(:, 28);
strain_30_30 = data_30_30(:, 25);  stress_30_30 = data_30_30(:, 28);
strain_0_60 = data_0_60(:, 25);   stress_0_60 = data_0_60(:, 28);

%% 2. VISUALIZATION
figure('Name','Hold Time Comparison','Color','w', 'Position', [100, 100, 900, 800]);
hold on;

% --- Add Zero-Axis Crosshairs ---
xline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off'); 
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');

% --- SHADED AREAS ---
fill(strain_60_0, stress_60_0, [0 0 0], 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(strain_30_30, stress_30_30, [1 0 0], 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(strain_0_60, stress_0_60, [0 0 1], 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% --- PLOT LINES & MARKERS ---
h1 = plot(strain_60_0, stress_60_0, '-s', 'Color', [0 0 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0 0 0], ...
    'MarkerIndices', round(linspace(1, length(strain_60_0), 100)));

h2 = plot(strain_30_30, stress_30_30, '-o', 'Color', [1 0 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], ...
    'MarkerIndices', round(linspace(1, length(strain_30_30), 100)));

h3 = plot(strain_0_60, stress_0_60, '-d', 'Color', [0 0 1], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0 0 1], ...
    'MarkerIndices', round(linspace(1, length(strain_0_60), 100)));

%% 3. ANNOTATION BOX (Strain Amplitude)
% Position is [x y width height] in normalized figure units
dim = [0.169 0.83 0.2 0.08]; 
str = {'\epsilon_a = 1.2%'};
annotation('textbox', dim, 'String', str, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 24, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

%% 4. FORMATTING
xlim([-0.02 0.02]); xticks([-0.02 -0.01 0 0.01 0.02]);
ylim([-1300 1300]); axis square;

xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontName', 'Times New Roman', ...
    'TickDir', 'in', 'Box', 'on');


%% 5. LEGEND
lgd = legend([h1, h2, h3], {'60/0 Hold', '30/30 Hold', '0/60 Hold'}, ...
    'Location', 'SouthEast', 'FontSize', 24,'FontName', 'Times New Roman');
set(lgd, 'FontName', 'Times New Roman', 'EdgeColor', 'k');

hold off;
fprintf('Plot generated with strain amplitude annotation box.\n');

% 2. Save the figure to a PDF
exportgraphics(gca, 'myPlot.pdf', 'ContentType', 'vector')