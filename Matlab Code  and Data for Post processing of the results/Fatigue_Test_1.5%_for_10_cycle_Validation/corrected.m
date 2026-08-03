clear all; close all; clc;

% =========================================================================
% 1. DATA LOADING
% =========================================================================
% Load Experimental data for Nickel Superalloy
if isfile('exp_1_5_.csv')
    data_exp = readmatrix('exp_1_5_.csv');
    strain_exp = data_exp(:, 1); 
    stress_exp = data_exp(:, 2);
else
    error('Experimental data file not found.');
end

% Load Simulated data (CPFEM / Pin Lu Model)
if isfile('Fatigue.csv')
    data_sim = readmatrix('Fatigue.csv');
    % Based on your MOOSE/CPFEM workflow
    strain_sim = data_sim(:, 25); % ε₃₃
    stress_sim = data_sim(:, 28); % σ₃₃
else
    warning('Simulated data (Fatigue.csv) not found.');
    strain_sim = []; stress_sim = [];
end

% =========================================================================
% 2. PLOT: HYSTERESIS LOOP COMPARISON
% =========================================================================
figure('Name','Fatigue Test: Simulated vs Experimental','Color','w');

% --- Plot Experimental data (Blue line with Upward Triangles) ---
if ~isempty(strain_exp)
    h1 = plot(strain_exp, stress_exp, '-^', ...
        'Color', 'b', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', 'b', ...
        'MarkerIndices', round(linspace(1, length(strain_exp), 25))); 
    hold on;
end

% --- Plot Simulated data (Magenta line with Downward Triangles) ---
if ~isempty(strain_sim)
    h2 = plot(strain_sim, stress_sim, '-v', ...
        'Color', [1 0 1], ... % Magenta
        'LineWidth', 1.8, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', [1 0 1], ...
        'MarkerIndices', round(linspace(1, length(strain_sim), 25))); 
end

% =========================================================================
% 3. PUBLICATION FORMATTING
% =========================================================================
% Typical limits for Superalloy Fatigue Cycles
xlim([-0.012 0.012]); 
ylim([-1200 1200]);

xlabel('Strain, \epsilon_{33}', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress, \sigma_{33} (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontName', 'Times New Roman', 'TickDir', 'in');
grid off; box on;

% --- Figure Label (a) ---
text(0.03, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 30, ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', 'VerticalAlignment', 'top');

% --- Legend (Exact Match to Reference) ---
% Using \epsilon_a to denote strain amplitude
legend([h1, h2], {'Exp. (\epsilon_a=1.0%)', 'Sim. (\epsilon_a=1.0%)'}, ...
    'Location', 'SouthEast', 'FontName', 'Times New Roman', 'FontSize', 22);
legend boxoff;

hold off;