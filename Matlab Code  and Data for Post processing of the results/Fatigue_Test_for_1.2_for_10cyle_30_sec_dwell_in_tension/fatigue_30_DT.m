clear all;
close all;
clc;

% --- Load Data ---
% Simulated data
data1 = readmatrix('Fatigue.csv');
strain_33_1 = data1(:, 25);   % ε₃₃ from first file
stress_33_1 = data1(:, 28);   % σ₃₃ from first file

% Experimental data
data2 = readmatrix('tension_dwell_1.2_.csv');
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
title('\bfFatigue Test: Simulated vs Experimental', 'FontSize', 16);
legend('\bfSimulated (Pin Lu Model)','\bfExperimental','Location','best');

grid off; box on;
set(gca, 'FontSize', 17, 'LineWidth', 1.4);
hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%                 NEW PLOTS ADDED BELOW                            %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Plot 4: New Plot 1 from Fatigue.csv ---
% !! SPECIFY YOUR COLUMNS HERE !!
new_x_col_1 = 1;  % <-- Change this to your desired X-axis column number
new_y_col_1 = 28;  % <-- Change this to your desired Y-axis column number
% Extract the new data
x_data_new_1 = data1(:, new_x_col_1);
y_data_new_1 = data1(:, new_y_col_1);
figure('Name','Stress vs Time','Color','w');
plot(x_data_new_1, y_data_new_1, '-', ...
    'LineWidth', 2, 'Color', 'r'); % Changed color to red
% !! UPDATE LABELS AND TITLE !!
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14, 'FontWeight', 'bold');
title('Stress vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% xlim([...]); % Optional: Set new limits if needed
% ylim([...]); % Optional: Set new limits if needed
% --- Plot 5: New Plot 2 from Fatigue.csv ---
% !! SPECIFY YOUR COLUMNS HERE !!
new_x_col_2 = 1;  % <-- Change this to your desired X-axis column number
new_y_col_2 = 22;  % <-- Change this to your desired Y-axis column number
% Extract the new data
x_data_new_2 = data1(:, new_x_col_2);
y_data_new_2 = data1(:, new_y_col_2);
figure('Name','Dislocation density vs Time','Color','w');
plot(x_data_new_2, y_data_new_2, '-', ...
    'LineWidth', 2, 'Color', [0.80 0.1740 0.9880]); % Changed color to green
% !! UPDATE LABELS AND TITLE !!
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Dislocation density', 'FontSize', 14, 'FontWeight', 'bold');
title('Dislocation density vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% xlim([...]); % Optional: Set new limits if needed
% ylim([...]); % Optional: Set new limits if needed
% plot6
% !! SPECIFY YOUR COLUMNS HERE !!
new_x_col_3 = 1;  % <-- Change this to your desired X-axis column number
new_y_col_3 = 21;  % <-- Change this to your desired Y-axis column number
% Extract the new data
x_data_new_3 = data1(:, new_x_col_3);
y_data_new_3 = data1(:, new_y_col_3);
figure('Name','Dislocation density vs Time','Color','w');
plot(x_data_new_3, y_data_new_3, '-', ...
    'LineWidth', 2, 'Color', [0.80 0.1740 0.9880]); % Changed color to green
% !! UPDATE LABELS AND TITLE !!
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('slip reistance', 'FontSize', 14, 'FontWeight', 'bold');
title('slip resistance vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
% xlim([...]); % Optional: Set new limits if needed
% ylim([...]); % Optional: Set new limits if needed



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%         NEW PLOT 6: COMBINED PLOT 4 and 5 (MODIFIED)             %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Stress and Dislocation Density vs Time', 'Color', 'w');

% --- MODIFIED: Define colors for clarity using user defaults ---
color_stress = [0 0.447 0.741];   % MATLAB blue
color_dislocation = [0.850 0.325 0.098]; % MATLAB orange-red

% --- Plot 1 (Stress) on Left Y-Axis ---
yyaxis left
plot(x_data_new_1, y_data_new_1, '-', 'LineWidth', 2, 'Color', color_stress);
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14, 'FontWeight', 'bold');
ax = gca;
ax.YColor = color_stress;

% --- Plot 2 (Dislocation) on Right Y-Axis ---
yyaxis right
plot(x_data_new_3, y_data_new_3, '-', 'LineWidth', 2, 'Color', color_dislocation);
ylabel('Slip resistance', 'FontSize', 14, 'FontWeight', 'bold');
ax.YColor = color_dislocation;

% --- Add Shared Labels and Title ---
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
title('Stress and Dislocation Density vs Time', 'FontSize', 16, 'FontWeight', 'bold');
legend('\bfStress (\sigma_{33})', '\bfDislocation Density', 'Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% --- MODIFIED: Added x-axis limit ---
xlim([0 675]);