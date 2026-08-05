clear all;
close all;
clc;

%%%%%%% Figure-2 (b)%%%%%%%%%%%%%%%%

% --- Load Data ---
% Simulated data
data1 = readmatrix('new_60_Tension.csv');
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

% --- Plot 1: Simulated (Original) ---
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

% --- Plot 2: Experimental (Original) ---
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

% =========================================================================
% =========================================================================
% PLOT 3: COMBINED (INCREASED BLACK MARKER DENSITY)
% =========================================================================
figure('Name','Fatigue Test: Simulated vs Experimental (High Marker Density)','Color','w');

% --- Simulated Data: Black Line + Squares (Increased Density) ---
% Increased 'linspace' to 50 for more markers
h1 = plot(strain_33_1, stress_33_1, '-s', ...
    'Color', [0 0 0], ...        % Black
    'LineWidth', 1.8, ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerIndices', round(linspace(1, length(strain_33_1), 150))); 
hold on;

% --- Experimental Data: Red Circles ---
h2 = plot(strain_33_2, stress_33_2, 'o', ...
    'MarkerEdgeColor', 'r', ...  
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 12, ...
    'LineWidth', 1.5, ...
    'MarkerIndices', round(linspace(1, length(strain_33_2), 200)));

% --- AXIS LIMITS & TICKS ---
xlim([-0.02 0.02]);
xticks([-0.02 -0.01 0 0.01 0.02]);
ylim([-1300 1300]); 
axis square;        

% --- LABELS ---
xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

% --- FORMATTING ---
set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontName', 'Times New Roman', ...
    'TickDir', 'in', 'FontWeight', 'normal');
grid off;
box on;

% --- Figure Label (b) ---
text(0.02, 0.98, '(b)', ...
    'Units', 'normalized', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');

% --- LEGEND ---
legend([h2, h1], {'Experimental', 'Simulated'}, ...
    'Location', 'SouthEast', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 24);
legend boxon;

hold off;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%                 NEW PLOTS ADDED BELOW                            %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Plot 4: New Plot 1 from Fatigue.csv ---
new_x_col_1 = 1;  
new_y_col_1 = 28;  
x_data_new_1 = data1(:, new_x_col_1);
y_data_new_1 = data1(:, new_y_col_1);

figure('Name','Stress vs Time','Color','w');
plot(x_data_new_1, y_data_new_1, '-', ...
    'LineWidth', 2, 'Color', 'r'); 
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\bf\sigma_{33} (Stress in MPa)', 'FontSize', 14, 'FontWeight', 'bold');
title('Stress vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);


% --- Plot 5: New Plot 2 from Fatigue.csv ---
new_x_col_2 = 1;  
new_y_col_2 = 22;  
x_data_new_2 = data1(:, new_x_col_2);
y_data_new_2 = data1(:, new_y_col_2);

figure('Name','Dislocation density vs Time','Color','w');
plot(x_data_new_2, y_data_new_2, '-', ...
    'LineWidth', 2, 'Color', [0.80 0.1740 0.9880]); 
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Dislocation density', 'FontSize', 14, 'FontWeight', 'bold');
title('Dislocation density vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);


% --- Plot 6: Slip Resistance ---
new_x_col_3 = 1; 
new_y_col_3 = 21;  
x_data_new_3 = data1(:, new_x_col_3);
y_data_new_3 = data1(:, new_y_col_3);

figure('Name','Slip resistance vs Time','Color','w');
plot(x_data_new_3, y_data_new_3, '-', ...
    'LineWidth', 2, 'Color', [0.80 0.1740 0.9880]); 
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Slip resistance', 'FontSize', 14, 'FontWeight', 'bold');
title('Slip resistance vs Time', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%         NEW PLOT 6: COMBINED PLOT 4 and 5 (MODIFIED)             %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Stress and Dislocation Density vs Time', 'Color', 'w');
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
xlabel('Time (sec)', 'FontSize', 14,'FontName', 'Times New Roman', 'FontWeight', 'bold');
title('Stress and Dislocation Density vs Time', 'FontSize', 16, 'FontWeight', 'bold');

% --- UPDATED LEGEND SETTINGS ---
legend({'\bfStress (\sigma_{33})', '\bfDislocation Density'}, ...
    'Location', 'SouthEast', ...       % Bottom-Right Corner
    'FontName', 'Times New Roman', ... % Font Change
    'FontSize', 12);

grid on; box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlim([0 675]);
% =========================================================================
% PLOT 3: COMBINED (WITH ZERO-AXIS GRID LINES)
% =========================================================================
figure('Name','Fatigue Test: Simulated vs Experimental (High Marker Density)','Color','w');
hold on;

% --- ADD ZERO-AXIS GRID LINES ---
% These draw the crosshairs at (0,0)
xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Vertical line at 0 strain
yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); % Horizontal line at 0 stress

% --- Simulated Data: Black Line + Squares (Increased Density) ---
h1 = plot(strain_33_1, stress_33_1, '-s', ...
    'Color', [0 0 0], ...        
    'LineWidth', 1.8, ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0], ...
    'MarkerIndices', round(linspace(1, length(strain_33_1), 150))); 

% --- Experimental Data: Red Circles ---
h2 = plot(strain_33_2, stress_33_2, 'o', ...
    'MarkerEdgeColor', 'b', ...  
    'MarkerFaceColor', 'b', ...
    'MarkerSize', 12, ...
    'LineWidth', 1.5, ...
    'MarkerIndices', round(linspace(1, length(strain_33_2), 200)));

% --- AXIS LIMITS & TICKS ---
xlim([-0.02 0.02]);
xticks([-0.02 -0.01 0 0.01 0.02]);
ylim([-1300 1300]); 
axis square;        

% --- LABELS ---
xlabel('Strain', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

% --- FORMATTING ---
set(gca, 'FontSize', 24, 'LineWidth', 2, 'FontName', 'Times New Roman', ...
    'TickDir', 'in', 'FontWeight', 'normal');
grid off;
box on;

% --- Figure Label (b) ---
text(0.02, 0.98, '(b)', ...
    'Units', 'normalized', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');

% --- LEGEND ---
legend([h2, h1], {'Experimental', 'Simulated'}, ...
    'Location', 'SouthEast', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 24);
legend boxon;
hold off;