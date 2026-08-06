close all;
clear all;
clc;

%%%%%%% Figure-6 (b)%%%%%%%%%%%%%5


%% ========================================================================
%% 1. INPUT DATA (Case-Specific Life Data)
%% ========================================================================
% Experimental Data (Cycles to Failure)
exp_00   = [66, 195, 450, 1400];       
exp_600  = [30, 73, 472, 713];         
exp_060  = [16, 176, 318, 534];        
exp_3030 = [10, 65, 220, 465];         

% Predicted Data (Simulation/Model)
pred_00   = [35, 111, 278, 1242];       
pred_600  = [17, 73, 472, 713];      
pred_060  = [9.4829, 273, 257, 736];       
pred_3030 = [6, 82, 311, 942.1];     

%% ========================================================================
%% 2. DEFINE CUSTOM COLORS (Publication Quality)
%% ========================================================================
c_grey   = [0.2, 0.2, 0.2];    % 0/0
c_red    = [0.85, 0.1, 0.1];   % 60/0
c_blue   = [0.0, 0.3, 0.7];    % 0/60
c_green  = [0.1, 0.6, 0.2];    % 30/30

%% ========================================================================
%% 3. SETUP PLOT LIMITS
%% ========================================================================
all_exp = [exp_00, exp_600, exp_060, exp_3030];
all_pred = [pred_00, pred_600, pred_060, pred_3030];

min_val = min([all_exp, all_pred]) * 0.1; 
max_val = max([all_exp, all_pred]) * 5;  
limits  = [min_val, max_val]; 

if isempty(limits) || isnan(limits(1)) || limits(1) <= 0
    limits = [1, 10000]; 
end

figure('Color', 'w'); hold on;

%% ========================================================================
%% 4. PLOT SCATTER BANDS
%% ========================================================================
% Perfect Agreement (1:1) - Black Solid
hPerf = plot(limits, limits, 'k-', 'LineWidth', 2); 

% +/- 2x Scatter - Dark Grey Dashed
h2x = plot(limits, limits * 2, 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
      plot(limits, limits / 2, 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

% +/- 3x Scatter - Muted Blue Dash-Dot
h3x = plot(limits, limits * 3, '-.', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.5);
      plot(limits, limits / 3, '-.', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.5);

%% ========================================================================
%% 5. PLOT DATA POINTS (Specific Markers & New Colors)
%% ========================================================================
% Case 0/0: Dark Grey Squares
h00 = loglog(exp_00, pred_00, 's', ...
    'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_grey);

% Case 60/0: Deep Red Circles
h600 = loglog(exp_600, pred_600, 'o', ...
    'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_red);

% Case 0/60: Royal Blue Triangles
h060 = loglog(exp_060, pred_060, '^', ...
    'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_blue);

% Case 30/30: Forest Green Diamonds
h3030 = loglog(exp_3030, pred_3030, 'd', ...
    'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_green);

%% ========================================================================
%% 6. FORMATTING (Times New Roman)
%% ========================================================================
axis square; grid off; box on;
xlim(limits); ylim(limits);
set(gca, 'XScale', 'log', 'YScale', 'log');

set(gca, 'LineWidth', 1.5, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal');

xlabel('Experimental Life (Cycles)', ...
    'FontSize', 26, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('Predicted Life (Cycles)', ...
    'FontSize', 26, 'FontName', 'Times New Roman', 'FontWeight', 'normal');

% --- LEGEND ---
lgd = legend([h00, h600, h060, h3030, hPerf, h2x, h3x], ...
    {'0/0', '60/0', '0/60', '30/30', ...
     'Perfect Match', '\pm 2x Scatter', '\pm 3x Scatter', '\pm 5x Scatter'}, ...
    'Location', 'northwest', 'NumColumns', 2);

lgd.FontSize = 24;
lgd.FontName = 'Times New Roman';
legend boxon;

% --- ANNOTATION ---
annotation('textbox', [0.65 0.15 0.25 0.12], ...
    'String', {'DD6, R_{\epsilon}=-1', 'Temperature = 760^{\circ}C'}, ...
    'Interpreter','tex', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth', 1.2);

hold off;
