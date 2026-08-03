clc;
clear;
close all;

%% 1. SETUP: DEFINE THE FIVE CASES
% =========================================================================
cases(1).filename = 'Fatigue_60_0_2_4.csv'; 
cases(1).epsilon_a = 0.012;      
cases(1).exp_life = 5;           

cases(2).filename = 'Fatigue_60_0_0_8_perc.csv'; 
cases(2).epsilon_a = 0.008;      
cases(2).exp_life = 600;         

cases(3).filename = 'Fatigue_60_0_1_8.csv'; 
cases(3).epsilon_a = 0.009;      
cases(3).exp_life = 297;          

cases(4).filename = 'Fatigue_60_0_2_0.csv'; 
cases(4).epsilon_a = 0.010;      
cases(4).exp_life = 100;         

cases(5).filename = 's_1_4.csv';  
cases(5).epsilon_a = 0.007;               
cases(5).exp_life = 6600;                 

%% GLOBAL MATERIAL CONSTANTS
n1 = 0.6;        
Q = 6.97e-19;    
T = 1033;        
R = 8.314;       
Sg = 0.055;      
Sc = 0.01268;   
Dfc = 0.055;  
q = 0.4; % Interaction exponent

% Initialize arrays
pred_lives_linear = zeros(1,5);
pred_lives_q      = zeros(1,5);
exp_lives         = zeros(1,5);

%% 2. MAIN COMPUTATION LOOP
for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:,1);
        entropy_rate_col = data(:,25);
        Delta_Sf = 0; Dc = 0;
        
        for i = 2:length(time_col)
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false; is_fatigue = false;
            
            % ----- CASE TIME WINDOWS -----
            if k==1
                if (t_current>552 && t_current<=612), is_creep=true;
                elseif (t_current>612 && t_current<=660), is_fatigue=true; end
            elseif k==2
                if (t_current>100 && t_current<=160), is_creep=true;
                elseif (t_current>160 && t_current<=200), is_fatigue=true; end
            elseif k==3
                if (t_current>297 && t_current<=357), is_creep=true;
                elseif (t_current>357 && t_current<=393), is_fatigue=true; end
            elseif k==4
                if (t_current>310 && t_current<=370), is_creep=true;
                elseif (t_current>370 && t_current<=410), is_fatigue=true; end
            elseif k==5
                if (t_current>7 && t_current<=37), is_creep=true;
                elseif (t_current>37 && t_current<=65), is_fatigue=true; end
            end
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + s_dot*dt;
            end
        end
        
        if Delta_Sf >= Sg
            Df = 1;
        else
            Df = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
        end
        
        % Linear Model: 1 / (Dc + Df)
        pred_lives_linear(k) = 1 / (Dc + Df);
        
        % Interaction Model: 1 / (Dc^q + Df^q)
        pred_lives_q(k) = 1 / (Dc^q + Df^q);
        
    else
        pred_lives_linear(k) = NaN;
        pred_lives_q(k)      = NaN;
    end
    exp_lives(k)  = cases(k).exp_life;
end

%% ========================================================================
%% 3. PLOT 1: LINEAR SUMMATION MODEL
%% ========================================================================
plotCorrelation(exp_lives, pred_lives_linear, 'Linear Summation Model', ...
    'N_f = (D_c + D_f)^{-1}');

%% ========================================================================
%% 4. PLOT 2: INTERACTION MODEL (q = 0.57)
%% ========================================================================
plotCorrelation(exp_lives, pred_lives_q, ['Interaction Model (q = ', num2str(q), ')'], ...
    ['N_f = (D_c^{', num2str(q), '} + D_f^{', num2str(q), '})^{-1}']);

%% ========================================================================
%% 5. PLOT 3: STRAIN RANGE vs LIFE (COMPARISON)
%% ========================================================================
strain_range = 2 * [cases.epsilon_a];
figure('Color','w','Position', [100 100 800 800]); hold on;

hExp = plot(exp_lives, strain_range, 'ok', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
hLin = plot(pred_lives_linear, strain_range, 'sr', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
hInt = plot(pred_lives_q, strain_range, 'db', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

set(gca, 'XScale', 'log', 'FontSize', 22, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlabel('Cycles to Failure, N_f', 'FontSize', 26);
ylabel('Strain Range, \Delta\epsilon', 'FontSize', 26);
grid on; box on; axis square;
legend([hExp, hLin, hInt], {'Experiment', 'Linear Model', 'Interaction Model'}, 'Location', 'northeast');
title('Model Comparison', 'FontWeight', 'normal');
hold off;

%% ========================================================================
%% HELPER FUNCTION FOR CORRELATION PLOTS
%% ========================================================================
function plotCorrelation(x, y, plotTitle, formulaStr)
    figure('Color','w'); hold on;
    colors  = {'r','b','g','m','k'}; 
    markers = {'s','d','^','v','p'};
    all_vals = [x, y];
    limits = [min(all_vals)*0.5, max(all_vals)*5];
    
    loglog(limits, limits, 'k-', 'LineWidth', 2); % Perfect match
    loglog(limits, limits*2, 'k--', 'LineWidth', 1.5); % 2x
    loglog(limits, limits*0.5, 'k--', 'LineWidth', 1.5);
    loglog(limits, limits*3, 'b-.', 'LineWidth', 1.5); % 3x
    loglog(limits, limits/3, 'b-.', 'LineWidth', 1.5);
    
    for k = 1:5
        loglog(x(k), y(k), markers{k}, 'MarkerSize', 12, 'LineWidth', 1.5, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
    end
    
    grid on; axis square; box on;
    xlim(limits); ylim(limits);
    set(gca, 'LineWidth', 1.5, 'FontSize', 22, 'FontName', 'Times New Roman');
    xlabel('Experimental Life (Cycles)');
    ylabel('Predicted Life (Cycles)');
    title(plotTitle, 'FontWeight', 'normal');
    
    annotation('textbox', [0.6 0.15 0.3 0.1], 'String', {formulaStr}, ...
        'Interpreter', 'tex', 'FontSize', 18, 'FontName', 'Times New Roman', ...
        'BackgroundColor', 'white');
end
%% 4. PLOT 2: STRAIN RANGE vs LIFE (WITH POLYNOMIAL FIT)
%% ========================================================================
strain_range = 2 * [cases.epsilon_a];
X_sim = pred_lives_linear; % Using the Interaction Model results
X_exp = exp_lives;
Y_vals = strain_range; 

figure('Color','w','Position', [900 100 800 800]); hold on;

% --- 1. POLYNOMIAL FIT FOR SIMULATION DATA ---
valid_idx = ~isnan(X_sim); 
logX = log10(X_sim(valid_idx));
Y_fit_data = Y_vals(valid_idx);

% Fit: log10(Nf) vs Strain Range
p = polyfit(logX, Y_fit_data, 2); 

% Generate points for a smooth dashed line
x_smooth = logspace(log10(min(X_sim(valid_idx))), log10(max(X_sim(valid_idx))), 100);
y_smooth = polyval(p, log10(x_smooth));

% Plot the Blue Dashed Fit Line
hFit = plot(x_smooth, y_smooth, 'b--', 'LineWidth', 2);

% --- 2. PLOT DATA POINTS ---
hExp = plot(X_exp, Y_vals, 'o', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
hSim = plot(X_sim, Y_vals, 's', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');

% --- FORMATTING ---
set(gca, 'XScale', 'log', 'FontSize', 22, 'LineWidth', 1.5, ...
         'FontName', 'Times New Roman', 'TickDir', 'in'); 
ylim([0.008 0.026]); % Adjusted to better fit your 0.007-0.012 epsilon_a range
xlim([1 100000]); 
xlabel('Cycles to Failure, N_f', 'FontSize', 26);
ylabel('Strain Range, \Delta\epsilon', 'FontSize', 26);
grid on; box on; axis square;
% --- ANNOTATION (Bottom Left) ---
annotation('textbox', [0.15 0.17 0.3 0.14], ...
    'String', {'DD6, R_{\epsilon}=-1, (60/0)', '■ 760°C'}, ...
    'Interpreter','tex', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth', 1.2);

% Label (c)
text(0.02, 0.98, '(c)', 'Units', 'normalized', 'FontSize', 24, ...
    'FontName', 'Times New Roman', 'VerticalAlignment', 'top');
% Legend and Annotation
lgd = legend([hExp, hSim, hFit], {'Experiment', 'Interaction Model', 'Sim. Trendline'}, ...
    'Location', 'northeast');
lgd.FontSize = 20;

hold off;