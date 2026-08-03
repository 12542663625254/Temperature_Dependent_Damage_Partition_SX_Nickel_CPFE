clear; clc; close all;

%% 1. GRID & DATA SETUP
% Define experimental parameters
eps_amps = [0.008, 0.009, 0.010, 0.012];
hold_times = [0, 60, 120, 240, 300];

% Synthetic damage ratio data (Z = Dc / (Dc + Df))
% Rows = Hold Times, Columns = Strain Amplitudes
Z_data_orig = [
    0.05,  0.05,  0.05,  0.05;  % Hold 0s   (Pure Fatigue)
    0.45,  0.45,  0.20,  0.80;  % Hold 60s  
    0.55,  0.30,  0.15,  0.85;  % Hold 120s 
    0.60,  0.50,  0.15,  0.90;  % Hold 240s
    0.80,  0.75,  0.30,  0.95   % Hold 300s (High Creep)
];

% Transpose the data matrix so X = Hold Time and Y = Strain Amplitude
Z_data = Z_data_orig';

% Create the base meshgrid
[X_mesh, Y_mesh] = meshgrid(hold_times, eps_amps);

%% 2. HIGH-FIDELITY INTERPOLATION
% Create a fine 500x500 grid for smooth, non-pixelated contours
[X_fine, Y_fine] = meshgrid(linspace(min(hold_times), max(hold_times), 500), ...
                            linspace(min(eps_amps), max(eps_amps), 500));

% Apply Modified Akima (makima) interpolation to prevent artificial oscillations
Z_fine = interp2(X_mesh, Y_mesh, Z_data, X_fine, Y_fine, 'makima');

%% 3. VISUALIZATION: CONTOURS & BOUNDARY
% Initialize a large, white-background figure
figure('Color', 'w', 'Position', [100 100 1200 700]);
ax = gca; hold on; box on;

% Shrink the main plot width to 60% to leave room for the custom legend
ax.Position = [0.10 0.12 0.55 0.80]; 

% A. Smooth Gradient Contour
contourf(X_fine, Y_fine, Z_fine, 150, 'LineStyle', 'none');

% B. Custom Diverging Colormap (Blue -> White -> Red)
cpts = [0 0 1; 1 1 1; 1 0 0]; 
cmap = interp1([0, 0.5, 1], cpts, linspace(0, 1, 256));
colormap(cmap); 
clim([0 1]); % Constrain color limits from 0 to 1

% C. Critical Mechanism Boundary Line (Dc = Df, or Z = 0.5)
contour(X_fine, Y_fine, Z_fine, [0.5, 0.5], 'k--', 'LineWidth', 3);

%% 4. REGIME ANNOTATIONS (TEXT BOXES)
% Coordinates for annotations: text(X, Y, 'String', ...)
% Southwest Box (Fatigue)
text(20, 0.0096, 'Fatigue Dominated', 'FontSize', 18, 'FontWeight', 'bold', ...
    'Color', [0 0 0.5], 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'Margin', 5, 'FontName', 'Times New Roman');

% Northeast Box (Creep)
text(260, 0.0115, 'Creep Dominated', 'FontSize', 18, 'FontWeight', 'bold', ...
    'Color', [0.5 0 0], 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'Margin', 5, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');

% Southeast Box (Creep)
text(280, 0.0085, 'Creep Dominated', 'FontSize', 18, 'FontWeight', 'bold', ...
    'Color', [0.5 0 0], 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'Margin', 5, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');

%% 5. FORMATTING, AXES, & LABELS
% Axis Labels & Limits
xlabel('Tensile Hold Time, $t_h$ (s)', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('Strain Amplitude, $\varepsilon_a$', 'Interpreter', 'latex', 'FontSize', 22);
xticks(hold_times);
yticks(eps_amps);
xlim([0, 300]); 
ylim([0.008, 0.012]);

% General Axis Properties
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on');

%% 6. COLORBAR & CUSTOM LEGEND
% Setup standard colorbar for reference
cb = colorbar;
cb.Position = [0.67 0.12 0.02 0.80]; % Place colorbar immediately to the right of the plot
cb.Ticks = [0, 0.25, 0.5, 0.75, 1.0];
cb.TickLabels = {'0 (Pure Fatigue)', '0.25', '0.50', '0.75', '1.0 (Pure Creep)'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;
ylabel(cb, 'Damage Ratio, $Z$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 270, 'VerticalAlignment', 'bottom');

% Create Custom Legend entries using off-screen patches
h_line    = plot(NaN, NaN, 'k--', 'LineWidth', 3);
h_fatigue = patch(NaN, NaN, [0.2 0.2 1], 'EdgeColor', 'k');      % Blue Swatch
h_mixed   = patch(NaN, NaN, [0.95 0.95 0.95], 'EdgeColor', 'k'); % Gray/White Swatch
h_creep   = patch(NaN, NaN, [1 0.2 0.2], 'EdgeColor', 'k');      % Red Swatch

% Compile legend
lgd = legend([h_line, h_fatigue, h_mixed, h_creep], ...
    {'Mechanism Boundary ($D_c = D_f$)', 'Fatigue Dominated ($D_c < D_f$)', ...
     'Mixed/Transitional ($D_c \approx D_f$)', 'Creep Dominated ($D_c > D_f$)'}, ...
     'Interpreter', 'latex', 'FontSize', 16, 'Location', 'eastoutside');

% Shift legend far to the right so it clears the colorbar completely
lgd.Position(1) = 0.73; 
lgd.Position(2) = 0.40; % Centers it vertically
lgd.Box = 'on';
lgd.LineWidth = 1.0;

hold off;