close all;
clear all;
clc;


%%%%%%%%%%%%% Figure-6(a)%%%%%%%%%%%%%%%%%



%% ========================================================================
%% 1. INPUT DATA (Case-Specific Life Data)
%% ========================================================================
% Experimental Data (Cycles to Failure)
exp_00   = [30, 285, 850, 3557];       
exp_600  = [5, 100, 297, 600];         
exp_060  = [16, 176, 318, 534];        
exp_3030 = [5, 250, 370, 562];         

% Predicted Data (Simulation/Model)
pred_00   = [59, 275, 685.1, 2381];       
pred_600  = [9.6467, 140.4, 291, 721.85];      
pred_060  = [9.4829, 273, 257, 736];       
pred_3030 = [5.0002, 130, 205, 524.44];     

%% ========================================================================
%% 2. DEFINE CUSTOM COLORS (Publication Quality)
%% ========================================================================
c_grey   = [0.2, 0.2, 0.2];    % 0/0
c_red    = [0.85, 0.1, 0.1];   % 60/0
c_blue   = [0.0, 0.3, 0.7];    % 0/60
c_green  = [0.1, 0.6, 0.2];    % 30/30

%% ========================================================================
%% FIGURE 1: EXPERIMENTAL vs PREDICTED LIFE (Correlation Plot)
%% ========================================================================
all_exp  = [exp_00, exp_600, exp_060, exp_3030];
all_pred = [pred_00, pred_600, pred_060, pred_3030];

min_val = min([all_exp, all_pred]) * 0.1; 
max_val = max([all_exp, all_pred]) * 5;  
limits  = [min_val, max_val]; 

if isempty(limits) || isnan(limits(1)) || limits(1) <= 0
    limits = [1, 10000]; 
end

figure('Color', 'w'); hold on;

% Perfect Agreement (1:1) - Black Solid
hPerf = plot(limits, limits, 'k-', 'LineWidth', 2); 

% +/- 2x Scatter
h2x = plot(limits, limits * 2, 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
      plot(limits, limits / 2, 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

% +/- 3x Scatter
h3x = plot(limits, limits * 3, '-.', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.5);
      plot(limits, limits / 3, '-.', 'Color', [0.3 0.5 0.9], 'LineWidth', 1.5);

% Case 0/0: Dark Grey Squares
h00 = loglog(exp_00, pred_00, 's', 'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_grey);

% Case 60/0: Deep Red Circles
h600 = loglog(exp_600, pred_600, 'o', 'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_red);

% Case 0/60: Royal Blue Triangles
h060 = loglog(exp_060, pred_060, '^', 'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_blue);

% Case 30/30: Forest Green Diamonds
h3030 = loglog(exp_3030, pred_3030, 'd', 'MarkerSize', 11, 'LineWidth', 1.2, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_green);

% --- FORMATTING ---
axis square; grid off; box on;
xlim(limits); ylim(limits);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'LineWidth', 1.5, 'FontSize', 24, ...
    'FontName', 'Times New Roman', 'FontWeight', 'normal');

xlabel('Experimental Life (Cycles)', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Predicted Life (Cycles)', 'FontSize', 26, 'FontName', 'Times New Roman');

lgd = legend([h00, h600, h060, h3030, hPerf, h2x, h3x], ...
    {'0/0', '60/0', '0/60', '30/30', ...
     'Perfect Match', '\pm 2x Scatter', '\pm 3x Scatter'}, ...
    'Location', 'northwest', 'NumColumns', 2);
lgd.FontSize = 24;
lgd.FontName = 'Times New Roman';
legend boxon;

annotation('textbox', [0.65 0.15 0.25 0.12], ...
    'String', {'DD6, R_{\epsilon}=-1', 'Temperature = 760^{\circ}C'}, ...
    'Interpreter','tex', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'BackgroundColor','white', 'EdgeColor','black', 'LineWidth', 1.2);

hold off;
