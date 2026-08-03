clear all; close all; clc;

% =========================================================================
% 1. LOAD DATA
% =========================================================================
% Load Experimental data
if isfile('fatigue_test_1_.csv')
    data_exp = readmatrix('fatigue_test_1_.csv');
    strain_exp = data_exp(:, 1); 
    stress_exp = data_exp(:, 2);
else
    error('Experimental data file not found.');
end

% Load Simulated data
if isfile('Fatigue.csv')
    data_sim = readmatrix('Fatigue.csv');
    strain_sim = data_sim(:, 25);
    stress_sim = data_sim(:, 28);
else
    warning('Simulated data (Fatigue.csv) not found.');
    strain_sim = []; stress_sim = [];
end

% =========================================================================
% 2. PLOT: COMBINED HYSTERESIS LOOP
% =========================================================================
figure('Name','Fatigue Test: Simulated vs Experimental','Color','w');

% --- Plot Experimental data (Blue line with Upward Triangles) ---
if ~isempty(strain_exp)
    h1 = plot(strain_exp, stress_exp, '-^', ...
        'Color', 'b', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', ...
        'MarkerIndices', round(linspace(1, length(strain_exp), 20))); % Sparse markers for clarity
    hold on;
end

% --- Plot Simulated data (Magenta line with Downward Triangles) ---
if ~isempty(strain_sim)
    h2 = plot(strain_sim, stress_sim, '-v', ...
        'Color', [1 0 1], ... % Magenta
        'LineWidth', 1.5, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', [1 0 1], ...
        'MarkerIndices', round(linspace(1, length(strain_sim), 20))); 
end

% --- Formatting ---
xlim([-0.012 0.012]); 
ylim([-1200 1200]);
xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

set(gca, 'FontSize', 24, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'TickDir', 'in');
grid off; box on;

% --- Figure Label (a) ---
text(0.03, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 28, ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal', 'VerticalAlignment', 'top');

% --- Legend (Formatted like your image) ---
% Note: \epsilon_a is used for strain amplitude as seen in your image
legend([h1, h2], {'Exp. (\epsilon_a=1.0%)', 'Sim. (\epsilon_a=1.0%)'}, ...
    'Location', 'SouthEast', 'FontName', 'Times New Roman', 'FontSize', 22);
legend boxoff;

hold off;