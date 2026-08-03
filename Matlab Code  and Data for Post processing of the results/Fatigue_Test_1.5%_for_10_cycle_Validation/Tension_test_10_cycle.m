clear all;
close all;
clc;

% =========================================================================
% 1. LOAD DATA
% =========================================================================
% Simulated data
if isfile('Fatigue.csv')
    data1 = readmatrix('Fatigue.csv');
    strain_33_1 = data1(:, 25);   % ε₃₃ from first file
    stress_33_1 = data1(:, 28);   % σ₃₃ from first file
else
    error('Fatigue.csv not found.');
end

% Experimental data
if isfile('exp_1_5_.csv')
    data2 = readmatrix('exp_1_5_.csv');
    strain_33_2 = data2(:, 1);    % ε₃₃ from second file
    stress_33_2 = data2(:, 2);    % σ₃₃ from second file
else
    % Fallback for safety
    if isfile('creep_fatigue_test_data.csv')
        data2 = readmatrix('creep_fatigue_test_data.csv');
        strain_33_2 = data2(:, 1);
        stress_33_2 = data2(:, 2);
    else
        warning('Experimental data file not found.');
        strain_33_2 = []; stress_33_2 = [];
    end
end

% =========================================================================
% 2. COLORS AND STYLES
% =========================================================================
color1 = [0 0.447 0.741];   % MATLAB blue
color2 = [0.850 0.325 0.098]; % MATLAB orange-red
marker1 = 'o';
marker2 = 's';

% =========================================================================
% 3. PLOT 1: SIMULATED (Standard View)
% =========================================================================
figure('Name','Tension Test (Simulated)','Color','w');
plot(strain_33_1, stress_33_1, '-', ...
    'LineWidth', 2, 'MarkerSize', 1, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color1, ...
    'Color', color1);
xlim([0 40]);
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 18, 'FontName', 'Times New Roman');
title('\bfTest 1: Tension Test (Simulated)', 'FontSize', 16, 'FontName', 'Times New Roman');
grid off; box on;
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'FontName', 'Times New Roman');

% =========================================================================
% 4. PLOT 2: EXPERIMENTAL (Standard View)
% =========================================================================
figure('Name','Tension Test (Experimental)','Color','w');
if ~isempty(strain_33_2)
    plot(strain_33_2, stress_33_2, ['--' marker2], ...
        'LineWidth', 1, 'MarkerSize', 5, ...
        'MarkerFaceColor', color2, ...
        'MarkerEdgeColor', 'k', ...
        'Color', color2);
end
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 18, 'FontName', 'Times New Roman');
title('\bfTest 2: Tension Test (Experimental)', 'FontSize', 16, 'FontName', 'Times New Roman');
grid on; box on;
set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'FontName', 'Times New Roman');

% =========================================================================
% =========================================================================
% 5. PLOT 3: COMBINED (NO LABELS)
% =========================================================================
figure('Name','Fatigue Test: Combined','Color','w');

color_sim = 'b';      % Blue Line
color_exp = 'r';      % Red Markers

% Plotting the Simulated data (Blue Line)
h1 = plot(strain_33_1, stress_33_1, '-', 'LineWidth', 2.5, 'Color', color_sim); 
hold on;

% Plotting the Experimental data (Red Markers, Size 12)
if ~isempty(strain_33_2)
    h2 = plot(strain_33_2, stress_33_2, 'o', ...
        'LineWidth', 1.5, ...       
        'MarkerSize', 9, ...       % Bigger Markers
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', color_exp); % Red Markers
end

% Formatting
xlim([-0.025 0.025]); 
xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 24, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'TickDir', 'in');
grid off; box on;

% --- Figure Label (a) ---
text(0.03, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 24, ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal', 'VerticalAlignment', 'top');

legend([h1, h2], {'Simulated', 'Experimental'}, 'Location', 'SouthEast', ...
    'FontName', 'Times New Roman', 'FontSize', 24);
legend boxon;

% --- Indicating Fatigue (Width of the Loop) ---
%line([-0.01, 0.01], [0, 0], 'Color', 'b', 'LineStyle', '--');
%text(0, 50, 'Fatigue (\Delta\epsilon_p)', 'HorizontalAlignment', 'center', 'FontSize', 14);

hold off;
% =========================================================================
% ERROR CALCULATION AT SPECIFIC STRAIN POINT (-0.015)
% =========================================================================
target_strain = -0.015;

if ~isempty(strain_33_1) && ~isempty(strain_33_2)
    % 1. Find the index in Experimental data closest to -0.015
    [~, idx_exp] = min(abs(strain_33_2 - target_strain));
    exp_stress_val = stress_33_2(idx_exp);
    exp_strain_val = strain_33_2(idx_exp); % Actual strain reached
    
    % 2. Find the index in Simulated data closest to -0.015
    [~, idx_sim] = min(abs(strain_33_1 - target_strain));
    sim_stress_val = stress_33_1(idx_sim);
    
    % 3. Calculate Percentage Error
    point_error_015 = abs(exp_stress_val - sim_stress_val) / abs(exp_stress_val) * 100;

    % 4. Display results
    fprintf('\n--- Analysis at Strain: %.4f ---\n', target_strain);
    fprintf('Experimental Stress: %.2f MPa [at strain %.5f]\n', exp_stress_val, exp_strain_val);
    fprintf('Simulated Stress:    %.2f MPa\n', sim_stress_val);
    fprintf('Percentage Error:    %.2f%%\n', point_error_015);
    fprintf('------------------------------------------\n');
end
