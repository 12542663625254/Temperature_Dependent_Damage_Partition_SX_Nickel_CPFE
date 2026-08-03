clc; clear; close all;
%% ========================================================================
%% 1. INPUT DATA
%% ========================================================================
% Predicted Data (Simulation/Model)
pred_00   = [59, 275, 685.1, 2381,11500];       
pred_600  = [9.6467, 140.4, 291,721.85,6840 ];      
pred_060  = [9.4829, 273, 257, 736,2500];       
pred_3030 = [5.0002, 130, 205, 524.44,1300];  




strain_range_vals = [0.024, 0.020, 0.018, 0.016, 0.014]; 

%% ========================================================================
%% 2. COLORS & DATA STRUCTURE
%% ========================================================================
c_black   = [0.0, 0.0, 0.0];    
c_blue    = [0.0, 0.0, 1.0];    
c_magenta = [1.0, 0.0, 1.0];    
c_orange  = [1.0, 0.5, 0.0];    
c_red     = [1.0, 0.0, 0.0];    

p_data(1).sim = pred_00;   p_data(1).clr = c_black;   p_data(1).mkr = 's'; p_data(1).lbl = '(0/0)';
p_data(2).sim = pred_600;  p_data(2).clr = c_blue;    p_data(2).mkr = 'o'; p_data(2).lbl = '(60/0)';
p_data(3).sim = pred_060;  p_data(3).clr = c_magenta; p_data(3).mkr = '^'; p_data(3).lbl = '(0/60)';
p_data(4).sim = pred_3030; p_data(4).clr = c_orange;  p_data(4).mkr = 'd'; p_data(4).lbl = '(30/30)';

%% ========================================================================
%% FIGURE 3: SIMULATION DATA WITH SCATTER BAND
%% ========================================================================
figure('Color', 'w', 'Name', 'Global Band Simulation', 'Position', [200 200 900 850]); hold on;

all_preds = [pred_00; pred_600; pred_060; pred_3030];
sim_leg_handles = []; 

for i = 1:length(p_data)
    X_sim = p_data(i).sim;
    clr = p_data(i).clr;
    mkr = p_data(i).mkr;
    
    p = polyfit(log10(X_sim), strain_range_vals, 2); 
    x_smooth = logspace(log10(0.1), log10(100000), 300);
    plot(x_smooth, polyval(p, log10(x_smooth)), '-', 'Color', clr, ...
        'LineWidth', 2.0, 'HandleVisibility', 'off');
    
    sim_leg_handles(i) = plot(X_sim, strain_range_vals, mkr, 'MarkerSize', 14, ...
        'LineWidth', 1.5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr);
end

%% ========================================================================
%% 3. DRAW THE GLOBAL BAND
%% ========================================================================
all_x_flat = all_preds(:);
all_y_flat = repmat(strain_range_vals, size(all_preds,1), 1); 
all_y_flat = all_y_flat(:);
log_x_midpoints = (log10(pred_00) + log10(pred_3030)) / 2;
p_mid = polyfit(log_x_midpoints, strain_range_vals, 1); 

rightward_nudge = 1.; 
p_mid(2) = p_mid(2) - p_mid(1) * log10(rightward_nudge);

x_center_vals = 10.^((all_y_flat - p_mid(2)) / p_mid(1));
req_factors = max(all_x_flat ./ x_center_vals, x_center_vals ./ all_x_flat).^2;
global_factor = ceil(max(req_factors) * 10) / 10; 

x_band = logspace(log10(0.1), log10(100000), 300);
plot(x_band, polyval(p_mid, log10(x_band * sqrt(global_factor))), 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
plot(x_band, polyval(p_mid, log10(x_band / sqrt(global_factor))), 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');

%% ========================================================================
%% 4. VISUALIZE FACTOR ARROW & TEXT
%% ========================================================================
y_viz = 0.019;
m = p_mid(1); c = p_mid(2);
log_x_mid = (y_viz - c) / m;
x_left_band = 10^(log_x_mid) / sqrt(global_factor);
x_right_band = 10^(log_x_mid) * sqrt(global_factor);

line([x_left_band, x_right_band], [y_viz, y_viz], ...
    'LineStyle', '--', 'Color', c_red, 'LineWidth', 2.5, 'HandleVisibility', 'off');

plot(x_left_band, y_viz, '<', 'MarkerEdgeColor', c_red, 'MarkerFaceColor', c_red, ...
    'MarkerSize', 10, 'HandleVisibility', 'off');
plot(x_right_band, y_viz, '>', 'MarkerEdgeColor', c_red, 'MarkerFaceColor', c_red, ...
    'MarkerSize', 10, 'HandleVisibility', 'off');

% Position text to the right, outside the band line
text(x_right_band * 1.1, y_viz + 0.0002, ...
    sprintf('%.1fx Scatter', global_factor), ...
    'FontSize', 22, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
    'Color', 'k', 'FontName', 'Times New Roman');

%% ========================================================================
%% 5. FORMATTING & EXPORT
%% ========================================================================
% Moved legend to southwest to avoid line interference
lgd_sim = legend(sim_leg_handles, {p_data.lbl}, 'Location', 'northeast');
lgd_sim.FontSize = 24; title(lgd_sim, 'Hold (s/s)'); 

xlim([1 50000]); ylim([0.014 0.024]);
set(gca, 'XScale', 'log', 'FontSize', 22, 'LineWidth', 2.0, ...
    'FontName', 'Times New Roman', 'TickDir', 'in'); 
box on; axis square; grid off;

xlabel('Cycles to Failure, $N_f$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Strain Range, $\Delta\epsilon$', 'Interpreter', 'latex', 'FontSize', 26);

% Updated annotation with 'scatter band' on the next line
annotation('textbox', [0.15 0.15 0.30 0.13], ...
    'String', {'Temp = 760^{\circ}C', ...
               sprintf('Red solid: %.1fx', global_factor), ...
               'scatter band'}, ...
    'FontSize', 20, 'FontName', 'Times New Roman', ...
    'BackgroundColor', 'white', 'EdgeColor', 'k', 'LineWidth', 1.0);
% Final export to PDF
exportgraphics(gcf, 'Scatter_Band_Plot.pdf', 'ContentType', 'vector');
hold off;