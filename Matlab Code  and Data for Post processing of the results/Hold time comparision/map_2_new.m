clear; clc; close all;

%% 1. COORDINATE SETUP (From your grid)
eps_amps_sorted = [0.008, 0.009, 0.010, 0.012];
hold_times_sorted = [0, 60, 120, 240, 300];

% Create Mesh
[X_mesh, Y_mesh] = meshgrid(hold_times_sorted, eps_amps_sorted);

%% 2. DATA PROCESSING (Replace with your actual all_Dc / all_Df)
% Z = Dc / (Dc + Df)
% Boundary logic: Dc = 0.5*Df -> Z = 0.333 | Df = 0.5*Dc -> Z = 0.667
Z_ratio = zeros(size(X_mesh)); 

% Placeholder loop for your variables:
% for d = 1:5, for k = 1:4 ... Z_ratio(e_idx, h_idx) = dc/(dc+df); end, end
% --- MOCK DATA FOR DEMONSTRATION ---
Z_ratio = [0.00, 0.45, 0.65, 0.85, 0.95;
           0.00, 0.30, 0.50, 0.70, 0.85;
           0.00, 0.15, 0.35, 0.55, 0.65;
           0.00, 0.05, 0.15, 0.30, 0.45];

%% 3. HIGH-FIDELITY INTERPOLATION
% Using 'makima' for smooth, non-oscillatory mechanism boundaries
[X_fine, Y_fine] = meshgrid(linspace(min(hold_times_sorted), max(hold_times_sorted), 500), ...
                            linspace(min(eps_amps_sorted), max(eps_amps_sorted), 500));
Z_fine = interp2(X_mesh, Y_mesh, Z_ratio, X_fine, Y_fine, 'makima');

%% 4. VISUALIZATION (ENHANCED INTENSITY)
figure('Color', 'w', 'Position', [100 100 1100 800]);
ax = gca; hold on; box on;

% Background: Smooth Diverging Contours
% Using 'shading interp' after contourf for maximum color density
[~, hF] = contourf(X_fine, Y_fine, Z_fine, 200, 'LineStyle', 'none');
shading interp; 

% Critical Regime Boundaries
[C, h] = contour(X_fine, Y_fine, Z_fine, [0.333, 0.667], 'k--', 'LineWidth', 3);
clabel(C, h, 'FontSize', 14, 'FontName', 'Times New Roman', 'LabelSpacing', 500);

% --- VIBRANT DIVERGING COLORMAP ---
% Reduced the 'white' spread to increase color saturation at the poles
cpts = [0.0 0.0 0.8;  % Vivid Blue (Fatigue)
        0.9 0.9 1.0;  % Very Light Blue/White transition
        1.0 0.9 0.9;  % Very Light Red transition
        0.8 0.0 0.0]; % Vivid Red (Creep)

% Interpolate the map to be denser
cmap = interp1([0, 0.45, 0.55, 1], cpts, linspace(0, 1, 256));
colormap(cmap);
clim([0 1]); 

%% 5. REGIME ANNOTATIONS
% Using white or high-contrast text for visibility against intense colors
text(40, 0.011, {'FATIGUE','DOMINATED'}, 'Color', 'w', 'FontSize', 18, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
text(150, 0.009, 'MIXED DAMAGE', 'Color', 'k', 'FontSize', 18, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
text(250, 0.0085, {'CREEP','DOMINATED'}, 'Color', 'w', 'FontSize', 18, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');

%% 6. FORMATTING & COLORBAR
xlabel('Tensile Hold Time, $t_h$ (s)', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Strain Amplitude, $\varepsilon_a$', 'Interpreter', 'latex', 'FontSize', 24);
cb = colorbar;
cb.Ticks = [0, 0.333, 0.5, 0.667, 1];
cb.TickLabels = {'$D_c=0$', '$D_c=0.5D_f$', '$D_c=D_f$', '$D_f=0.5D_c$', '$D_f=0$'};
cb.TickLabelInterpreter = 'latex';
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 2);
xticks(hold_times_sorted);
yticks(eps_amps_sorted);
xlim([0, 300]);
ylim([0.008, 0.012]); 

grid off;

%% 7. SAVE FIGURE
% Use exportgraphics (Recommended for MATLAB R2020a and newer)
exportgraphics(gcf, 'Damage_Map_760C.pdf', 'ContentType', 'vector');