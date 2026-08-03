clear all; close all; clc;

% =========================================================================
% 1. POSITIONING LOGIC: SQUARE FIGURE DISPLAY
% =========================================================================
% Get current screen size
screen_size = get(0, 'Screensize'); 
% Define a square dimension (matches the size of your reference plot)
fig_dim = 750; 
% Center the square figure on the screen
left_pos = (screen_size(3) - fig_dim) / 2;
bottom_pos = (screen_size(4) - fig_dim) / 2;
square_pos = [left_pos, bottom_pos, fig_dim, fig_dim];

% =========================================================================
% 2. LOAD DATA
% =========================================================================
% Data Set 1 (1.0% Strain Amplitude)
if isfile('Fatigue.csv') && isfile('fatigue_test_1_.csv')
    data_sim_10 = readmatrix('Fatigue.csv'); 
    data_exp_10 = readmatrix('fatigue_test_1_.csv');
else
    error('Data files for 1.0%% strain amplitude not found.');
end

% Data Set 2 (1.5% Strain Amplitude)
if isfile('Fatigue_1_5.csv') && isfile('exp_1_5_.csv')
    data_sim_15 = readmatrix('Fatigue_1_5.csv'); 
    data_exp_15 = readmatrix('exp_1_5_.csv');
else
    error('Data files for 1.5%% strain amplitude not found.');
end

% =========================================================================
% 3. PLOT: COMBINED HYSTERESIS LOOPS
% =========================================================================
figure('Name','Creep-Fatigue Square Plot','Color','w','Position', square_pos);
hold on;

% --- 0. ADD ZERO-AXIS GRID LINES ---
% This draws the crosshairs at (0,0)
xline(0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5); % Vertical line at 0 strain
yline(0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.5); % Horizontal line at 0 stress

% --- 1.0% Strain Amplitude Data ---
% Exp: Blue Upward Triangles
h1 = plot(data_exp_10(:,1), data_exp_10(:,2), '^', ...
    'MarkerEdgeColor', 'b', 'MarkerSize', 12, 'MarkerFaceColor', 'b', ...
    'MarkerIndices', round(linspace(1, length(data_exp_10), 25)));

% Sim: Magenta Downward Triangles (Line + Markers)
h2 = plot(data_sim_10(:,25), data_sim_10(:,28), '-v', 'Color', [1 0 1], ...
    'LineWidth', 1.8, 'MarkerSize', 10, 'MarkerFaceColor', [1 0 1], ...
    'MarkerIndices', round(linspace(1, length(data_sim_10), 20)));

% --- 1.5% Strain Amplitude Data ---
% Exp: Red Circles
h3 = plot(data_exp_15(:,1), data_exp_15(:,2), 'o', ...
    'MarkerEdgeColor', 'r', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'MarkerIndices', round(linspace(1, length(data_exp_15), 25)));

% Sim: Black Squares (Line + Markers)
h4 = plot(data_sim_15(:,25), data_sim_15(:,28), '-s', 'Color', [0 0 0], ...
    'LineWidth', 1.8, 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0], ...
    'MarkerIndices', round(linspace(1, length(data_sim_15), 20)));

% =========================================================================
% 4. FORMATTING (SQUARE & PUBLICATION READY)
% =========================================================================
xlim([-0.02 0.02]); 
ylim([-1300 1300]);
axis square; % Critical for maintaining physical loop proportions

xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontName', 'Times New Roman', ...
         'TickDir', 'in', 'FontWeight', 'normal');

grid off; box on;

% --- Figure Label (a) Top-Left ---
text(0.02, 0.98, '(a)', ...
    'Units', 'normalized', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');

% --- Legend ---
lgd = legend([h1, h2, h3, h4], ...
    {'Exp. (\epsilon_a=1.0%)', 'Sim. (\epsilon_a=1.0%)', ...
     'Exp. (\epsilon_a=1.5%)', 'Sim. (\epsilon_a=1.5%)'}, ...
    'Location', 'SouthEast', 'FontName', 'Times New Roman', 'FontSize', 24);
legend boxon;

hold off;