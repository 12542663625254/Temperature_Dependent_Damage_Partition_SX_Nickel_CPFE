clear all;
close all;
clc;
% --- Load Data ---
% Original datasets
data1 = readmatrix('new.csv');
strain_33_1 = data1(:, 1);   % ε₃₃ from Sim 1
stress_33_1 = data1(:, 2);   % σ₃₃ from Sim 1
data2 = readmatrix('stress_exp_amp.csv');
strain_33_2 = data2(:, 1);    % ε₃₃ from Exp 1
stress_33_2 = data2(:, 2);    % σ₃₃ from Exp 1
% --- NEW: Load Second Set of Data ---
% !!! PLEASE UPDATE THESE FILENAMES to your new files !!!
data3 = readmatrix('1%_simulated_stress.csv'); 
strain_33_3 = data3(:, 1);   % ε₃₃ from Sim 2
stress_33_3 = data3(:, 2);   % σ₃₃ from Sim 2
data4 = readmatrix('stres_1%.csv');
strain_33_4 = data4(:, 1);    % ε₃₃ from Exp 2
stress_33_4 = data4(:, 2);    % σ₃₃ from Exp 2
% --- Colors and Styles ---
color1 = [0 0.447 0.741];   % MATLAB blue (Sim 1)
color2 = [0.850 0.325 0.098]; % MATLAB orange-red (Exp 1)
marker1 = 'o';
marker2 = 's';
% --- NEW: Colors and Styles for Second Set ---
color3 = [0.466 0.674 0.188]; % MATLAB green (Sim 2)
color4 = [0.494 0.184 0.556]; % MATLAB purple (Exp 2)
marker3 = '^'; % Triangle for Sim 2
marker4 = 'd'; % Diamond for Exp 2

% --- Plot 1: All Simulated ---
% This plot is unchanged
figure('Name','Tension Test (Simulated)','Color','w');
plot(strain_33_1, stress_33_1, '-', ...
    'LineWidth', 2, ...
    'Color', color1);
hold on;
plot(strain_33_3, stress_33_3, '--', ...
    'LineWidth', 2, ...
    'Color', color3);
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 14);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14);
title('\bfAll Simulated Tension Tests', 'FontSize', 16);
legend('\bfSimulated 1', '\bfSimulated 2', 'Location', 'best');
grid off; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold off;

% --- Plot 2: All Experimental ---
% This plot is unchanged
figure('Name','Tension Test (Experimental)','Color','w');
plot(strain_33_2, stress_33_2, ['--' marker2], ...
    'LineWidth', 1, 'MarkerSize', 5, ...
    'MarkerFaceColor', color2, ...
    'MarkerEdgeColor', 'k', ...
    'Color', color2);
hold on;
plot(strain_33_4, stress_33_4, [':' marker4], ...
    'LineWidth', 1, 'MarkerSize', 5, ...
    'MarkerFaceColor', color4, ...
    'MarkerEdgeColor', 'k', ...
    'Color', color4);
xlabel('\bf\epsilon_{33} (Strain)', 'FontSize', 14);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14);
title('\bfAll Experimental Tension Tests', 'FontSize', 16);
legend('\bfExperimental 1', '\bfExperimental 2', 'Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold off;

% --- Plot 3: Combined (with your new legend) ---
% This section has been updated to match your screenshot
figure('Name','Fatigue Test: Simulated vs Experimental','Color','w');

% --- MODIFIED ---
% Exp. (εa=1.0%) - Mapped to data2 (strain_33_2, stress_33_2)
% Red square, solid line (was black)
plot(strain_33_2, stress_33_2, '-sr', ...
    'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'DisplayName', 'Exp. (\epsilon_a=1.2%)');
hold on;

% Sim. (εa=1.0%) - Mapped to data1 (strain_33_1, stress_33_1)
% Red circle, solid line (unchanged)
plot(strain_33_1, stress_33_1, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'r', ...
    'Color', 'r', ...
    'DisplayName', 'Sim. (\epsilon_a=1.2%)');

% Exp. (εa=1.2%) - Mapped to data4 (strain_33_4, stress_33_4)
% Blue triangle, solid line (unchanged)
plot(strain_33_4, stress_33_4, '-^', ...
    'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'b', ...
    'Color', 'b', ...
    'DisplayName', 'Exp. (\epsilon_a=1.0%)');

% --- MODIFIED ---
% Sim. (εa=1.2%) - Mapped to data3 (strain_33_3, stress_33_3)
% Blue inverted triangle, solid line (was green)
plot(strain_33_3, stress_33_3, '-v', ...
    'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'b', ...
    'Color', 'b', ...
    'DisplayName', 'Sim. (\epsilon_a=1.0%)');

% Labels and title (using your original limits)
xlim([0 12]);
ylim([750 1200]);
xlabel('\bf Number of cycles', 'FontSize', 14);
ylabel('\bf Stress amplitude (MPa)', 'FontSize', 14);
title('\bfTension Test: Simulated vs Experimental', 'FontSize', 16);

% --- NEW LEGEND COMMAND ---
% This automatically uses the 'DisplayName' from each plot
legend('Location','northwest', 'Interpreter', 'tex', 'FontSize', 12, 'box', 'on');

grid off; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.4);
hold off;
% --- Plot 3: Combined Validation (1.2% Close Alignment) ---
figure('Name','Fatigue Test: Validation at 1.2%','Color','w','Position', [100 100 800 650]);

% 1. Experimental (εa=1.2%) - Red Squares
% Using a slightly larger marker to serve as the background 'anchor'
plot(strain_33_2, stress_33_2, 'sr', ...
    'LineWidth', 1.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', ... 
    'DisplayName', 'Exp. (\epsilon_a=1.2%)');
hold on;

% 2. Simulated (εa=1.2%) - Red Circles
% Plotted with a solid line to show the continuous predictive path
plot(strain_33_1, stress_33_1, '-or', ...
    'LineWidth', 2.0, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'w', ... 
    'MarkerEdgeColor', 'r', ...
    'DisplayName', 'Sim. (\epsilon_a=1.2%)');

% 3. Exp. (εa=1.0%) - Blue Triangles
plot(strain_33_4, stress_33_4, '^b', ...
    'LineWidth', 1.5, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Exp. (\epsilon_a=1.0%)');

% 4. Sim. (εa=1.0%) - Blue Inverted Triangles
plot(strain_33_3, stress_33_3, '-vb', ...
    'LineWidth', 2.0, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'b', ...
    'DisplayName', 'Sim. (\epsilon_a=1.0%)');

% --- Formatting for Publication ---
xlim([0 12]);
ylim([750 1200]);
xlabel('\bf Number of cycles', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('\bf Stress amplitude (MPa)', 'FontSize', 18, 'FontName', 'Times New Roman');

% Legend with white background for clarity
legend('Location','southeast', 'Interpreter', 'tex', 'FontSize', 14, ...
       'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', 'k');

grid on; box on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times New Roman', 'TickDir', 'in');

% Adding a text annotation for accuracy
text(0.05, 0.9, 'Error < 2.3%', 'Units', 'normalized', 'FontSize', 16, ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);

hold off;
% --- Plot 3: Stress Amplitude Evolution (Validation at 1.2% and 1.0%) ---
figure('Name','Figure 3: Stress Amplitude Validation','Color','w','Position', [100 100 1000 850]);

% 1. Experimental (εa=1.2%) - Red Squares
plot(strain_33_2, stress_33_2, 'sr', ...
    'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', ... 
    'DisplayName', 'Exp. (\epsilon_a=1.2%)');
hold on;

% 2. Simulated (εa=1.2%) - Red Circles (Solid Line)
plot(strain_33_1, stress_33_1, '-or', ...
    'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ... 
    'MarkerEdgeColor', 'r', ...
    'DisplayName', 'Sim. (\epsilon_a=1.2%)');

% 3. Experimental (εa=1.0%) - Blue Triangles
plot(strain_33_4, stress_33_4, '^b', ...
    'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Exp. (\epsilon_a=1.0%)');

% 4. Simulated (εa=1.0%) - Blue Inverted Triangles (Solid Line)
plot(strain_33_3, stress_33_3, '-vb', ...
    'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'b', ...
    'DisplayName', 'Sim. (\epsilon_a=1.0%)');

% --- PUBLICATION FORMATTING (Times New Roman, Large Fonts) ---
xlim([0 12]);
ylim([750 1200]);

% Labels with Font Size 26
xlabel('Number of cycles', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress amplitude (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

% Legend with Font Size 24
lgd = legend('Location', 'southeast', 'Interpreter', 'tex', 'FontSize', 24, ...
             'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', 'k');

% Axes Formatting
grid off; box on;
set(gca, 'FontSize', 24, ...                 % Set axis tick labels to 26
         'LineWidth', 2.0, ...              % Thicker box lines for clarity
         'FontName', 'Times New Roman', ... % Force Times New Roman
         'TickDir', 'in');

% Text Annotation for Accuracy
text(0.05, 0.92, 'Max Error < 1.7%', 'Units', 'normalized', ...
    'FontSize', 26, 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);
axis square;
hold off;

% --- Figure 3: Stress Amplitude Validation (Joined Experimental Lines) ---
figure('Name','Figure 3: Stress Amplitude Validation','Color','w','Position', [100 100 1000 850]);

% 1. Experimental (εa=1.2%) - Red Squares joined by Dashed Line
plot(strain_33_2, stress_33_2, '--sr', ...
    'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', ... 
    'DisplayName', 'Exp. (\epsilon_a=1.2%)');
hold on;

% 2. Simulated (εa=1.2%) - Red Circles joined by Solid Line
plot(strain_33_1, stress_33_1, '-or', ...
    'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ... 
    'MarkerEdgeColor', 'r', ...
    'DisplayName', 'Sim. (\epsilon_a=1.2%)');

% 3. Experimental (εa=1.0%) - Blue Triangles joined by Dashed Line
plot(strain_33_4, stress_33_4, '--^b', ...
    'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'k', ... 
    'DisplayName', 'Exp. (\epsilon_a=1.0%)');

% 4. Simulated (εa=1.0%) - Blue Inverted Triangles joined by Solid Line
plot(strain_33_3, stress_33_3, '-vb', ...
    'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', ... 
    'MarkerEdgeColor', 'b', ...
    'DisplayName', 'Sim. (\epsilon_a=1.0%)');

% --- PUBLICATION FORMATTING (Times New Roman, Large Fonts) ---
xlim([0 12]);
ylim([750 1200]);

% Labels with Font Size 26
xlabel('Number of cycles', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress amplitude (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

% Legend with Font Size 24
lgd = legend('Location', 'southeast', 'Interpreter', 'tex', 'FontSize', 24, ...
             'FontName', 'Times New Roman', 'Box', 'on', 'EdgeColor', 'k');

% Axes Formatting
grid off; box on;
set(gca, 'FontSize', 22, ...                 
         'LineWidth', 2.0, ...              
         'FontName', 'Times New Roman', ... 
         'TickDir', 'in');

% Text Annotation for Accuracy
text(0.05, 0.92, 'Max Error < 2.3%', 'Units', 'normalized', ...
    'FontSize', 24, 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);
axis square;
hold off;
% --- CALCULATE ERRORS ---

% 1. Error for epsilon_a = 1.2%
% Interpolate simulation (data1) onto experimental cycle points (data2)
sim_interp_12 = interp1(strain_33_1, stress_33_1, strain_33_2, 'linear', 'extrap');
error_12 = abs(stress_33_2 - sim_interp_12) ./ stress_33_2 * 100;
mean_error_12 = mean(error_12);
max_error_12 = max(error_12);

% 2. Error for epsilon_a = 1.0%
% Interpolate simulation (data3) onto experimental cycle points (data4)
sim_interp_10 = interp1(strain_33_3, stress_33_3, strain_33_4, 'linear', 'extrap');
error_10 = abs(stress_33_4 - sim_interp_10) ./ stress_33_4 * 100;
mean_error_10 = mean(error_10);
max_error_10 = max(error_10);

% Display results in Command Window
fprintf('--- Error Analysis ---\n');
fprintf('1.2%% Strain: Mean Error = %.2f%%, Max Error = %.2f%%\n', mean_error_12, max_error_12);
fprintf('1.0%% Strain: Mean Error = %.2f%%, Max Error = %.2f%%\n', mean_error_10, max_error_10);

% --- PLOT WITH CALCULATED ERROR ANNOTATION ---
figure('Name','Stress Amplitude Validation','Color','w','Position', [100 100 1000 850]);
hold on;

% 1. Experimental 1.2% (Dashed line)
plot(strain_33_2, stress_33_2, '--sr', 'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Exp. (\epsilon_a=1.2%)');

% 2. Simulated 1.2% (Solid line)
plot(strain_33_1, stress_33_1, '-or', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'DisplayName', 'Sim. (\epsilon_a=1.2%)');

% 3. Experimental 1.0% (Dashed line)
plot(strain_33_4, stress_33_4, '--^b', 'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Exp. (\epsilon_a=1.0%)');

% 4. Simulated 1.0% (Solid line)
plot(strain_33_3, stress_33_3, '-vb', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'DisplayName', 'Sim. (\epsilon_a=1.0%)');

% Publication Formatting
xlim([0 12]); ylim([750 1200]); axis square;
xlabel('Number of cycles', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress amplitude (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');
legend('Location', 'southeast', 'FontSize', 22, 'FontName', 'Times New Roman', 'Box', 'on');

% Dynamic Annotation using calculated Max Error
total_max_err = max(max_error_12, max_error_10);
text(0.05, 0.92, sprintf('Max Error: %.1f%%', total_max_err), 'Units', 'normalized', ...
    'FontSize', 24, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'Color', [0.2 0.5 0.2]);

set(gca, 'FontSize', 22, 'LineWidth', 2.0, 'FontName', 'Times New Roman', 'TickDir', 'in');
hold off;