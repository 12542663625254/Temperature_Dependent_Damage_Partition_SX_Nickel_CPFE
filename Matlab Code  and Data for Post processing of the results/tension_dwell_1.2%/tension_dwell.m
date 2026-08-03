clear all;
close all;
clc;

% --- Load Data ---
% Simulated data
data1 = readmatrix('Fatigue.csv');
strain_33_1 = data1(:, 25);   % ε₃₃ from first file
stress_33_1 = data1(:, 28);   % σ₃₃ from first file

% Experimental data
data2 = readmatrix('Dwell_Tension.csv');
strain_33_2 = data2(:, 1);    % ε₃₃ from second file
stress_33_2 = data2(:, 2);    % σ₃₃ from second file

% --- Colors and Styles ---
color1 = [0 0.447 0.741];   % MATLAB blue
color2 = [0.850 0.325 0.098]; % MATLAB orange-red
marker1 = 'o';
marker2 = 's';

% --- Plot 1: Simulated ---
figure('Name','Tension Test (Simulated)','Color','w');
plot(strain_33_1, stress_33_1, '-', ...
    'LineWidth', 2, 'MarkerSize', 1, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color1, ...
    'Color', color1);
ylim([2e8 200002000]);
xlim([0 40])
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 14);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14);
title('\bfTest 1: Tension Test (Simulated)', 'FontSize', 16);
grid off; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% --- Plot 2: Experimental ---
figure('Name','Tension Test (Experimental)','Color','w');
plot(strain_33_2, stress_33_2, ['--' marker2], ...
    'LineWidth', 1, 'MarkerSize', 5, ...
    'MarkerFaceColor', color2, ...
    'MarkerEdgeColor', 'k', ...
    'Color', color2);
ylim([2e8 2000002000]);
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 14);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14);
title('\bfTest 2: Tension Test (Experimental)', 'FontSize', 16);
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% --- Plot 3: Combined with new marker for simulated ---
figure('Name','Fatigue Test: Simulated vs Experimental','Color','w');

% Simulated data with triangle marker
plot(strain_33_1, stress_33_1, '-o', ...
    'LineWidth', 2, 'MarkerSize', 1, ...
    'MarkerFaceColor', color1, ...
    'MarkerEdgeColor', color1, ...
    'Color', color1);
hold on;

% Experimental data (keep same)
plot(strain_33_2, stress_33_2, 's', ...
    'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerFaceColor', color2, ...
    'MarkerEdgeColor', 'k', ...
    'Color', color2);

% Labels and title
xlim([-0.02 0.02]);
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 14);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14);
title('\bfTension Test: Simulated vs Experimental', 'FontSize', 16);
legend('\bfSimulated (Pin Lu Model)','\bfExperimental','Location','best');

grid off; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.4);
hold off;
