clc;
clear;
close all;

%% 1. SETUP: DEFINE THE FOUR CASES
% =========================================================================
% CASE 1: 1.2% Strain 
cases(1).filename = 'cyclic_2_4.csv'; 
cases(1).epsilon_a = 0.012;      
cases(1).exp_life = 30;           

% CASE 2: 0.8% Strain 
cases(2).filename = 'cyclic_2_0.csv'; 
cases(2).epsilon_a = 0.010;      
cases(2).exp_life = 290;         

% CASE 3: 0.9% Strain 
cases(3).filename = 'cyclic_1_8.csv'; 
cases(3).epsilon_a = 0.009;      
cases(3).exp_life = 850;          

% CASE 4: 1.0% Strain 
cases(4).filename = 'cyclic_1_6.csv'; 
cases(4).epsilon_a = 0.008;      
cases(4).exp_life = 3557;         

% (Case 5 Removed)

% =========================================================================
% GLOBAL MATERIAL CONSTANTS
Sg = 0.055;      
%Sc = 0.9 * Sg;
Sc = 0.01268;
D0 = 0;          
%Dfc = 0.9;
Dfc =0.04;
n1 = 0.6;        
Q = 6.97e-19;    
T = 1033;        
R = 8.314;       

% Pre-allocate for 4 cases
pred_lives = zeros(1, 4);
exp_lives = zeros(1, 4);

%% 2. MAIN PROCESSING LOOP
for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:, 1);          
        entropy_rate_col = data(:, 25); 
        
        Delta_Sf = 0; 
        Dc = 0;       
        num_points = length(time_col);
        
        for i = 2:num_points
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false;
            is_fatigue = false;
            
            % =========================================================
            % TIME LOGIC (Assuming Pure Fatigue for this script)
            % =========================================================
            if k == 1 
                if (t_current > 12 && t_current <= 60), is_fatigue = true; end
            elseif k == 2
               if (t_current > 10 && t_current <= 50), is_fatigue = true; end
            elseif k == 3
               if (t_current > 9 && t_current <= 45), is_fatigue = true; end
            elseif k == 4
                if (t_current > 8 && t_current <= 40), is_fatigue = true; end
            end
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1 - n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + (s_dot * dt);
            end
        end
        
        if Delta_Sf >= Sg
            Df = 1.0; 
        else
            numerator = Dfc - D0;
            denominator = log(1 - Sc/Sg);
            log_term = log(1 - Delta_Sf/Sg);
            Df = D0 + (numerator / denominator) * log_term;
        end
        
        pred_lives(k) = 1 / (Dc + Df);
    else
        pred_lives(k) = NaN;
    end
    exp_lives(k) = cases(k).exp_life;
end

%% ========================================================================
%% 3. PLOT 1: PREDICTED vs EXPERIMENTAL LIFE
%% ========================================================================
figure('Color', 'w');
colors = {'r', 'b', 'g', 'm'}; % Removed 'k' (black)
pred_markers = {'s', 'd', '^', 'v'}; % Removed '>'
exp_markers = {'x', '+', '*', 'p'}; % Removed 's'
all_vals = [exp_lives, pred_lives];
limits = [min(all_vals(all_vals>0)) * 0.1, max(all_vals) * 5];
if isempty(limits) || isnan(limits(1)), limits = [1 10000]; end

% --- Reference Lines ---
hPerf = loglog(limits, limits, 'k-', 'LineWidth', 2); hold on;
h2x = loglog(limits, limits * 2, 'k--', 'LineWidth', 1.5);
loglog(limits, limits * 0.5, 'k--', 'LineWidth', 1.5);
h3x = loglog(limits, limits * 3, 'b-.', 'LineWidth', 1.5);
loglog(limits, limits / 3, 'b-.', 'LineWidth', 1.5);

% --- Plot Data ---
for k = 1:4
    hPred(k) = loglog(exp_lives(k), pred_lives(k), pred_markers{k}, ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
    
    hExp(k) = loglog(exp_lives(k), exp_lives(k), exp_markers{k},...
        'MarkerSize',14,'LineWidth',2,'MarkerEdgeColor',colors{k},'MarkerFaceColor','none');
end

grid off; axis square;
xlim(limits); ylim(limits);

% --- Font Settings ---
set(gca, 'LineWidth', 1.5, 'FontSize', 22, 'FontName', 'Times New Roman');
xlabel('Experimental Life (Cycles)', 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');
ylabel('Predicted Life (Cycles)', 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');
title('Creep-Fatigue Life Prediction', 'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');

% --- Legend Construction ---
e1 = plot(nan,nan,exp_markers{1},'MarkerSize',11,'MarkerEdgeColor','r','LineWidth',1.8,'MarkerFaceColor','none');
e2 = plot(nan,nan,exp_markers{2},'MarkerSize',11,'MarkerEdgeColor','b','LineWidth',1.8,'MarkerFaceColor','none');
e3 = plot(nan,nan,exp_markers{3},'MarkerSize',11,'MarkerEdgeColor','g','LineWidth',1.8,'MarkerFaceColor','none');
e4 = plot(nan,nan,exp_markers{4},'MarkerSize',11,'MarkerEdgeColor','m','LineWidth',1.8,'MarkerFaceColor','none');

p1 = plot(nan,nan,pred_markers{1},'MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k');
p2 = plot(nan,nan,pred_markers{2},'MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k');
p3 = plot(nan,nan,pred_markers{3},'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');
p4 = plot(nan,nan,pred_markers{4},'MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k');

lgd = legend([e1 e2 e3 e4 p1 p2 p3 p4], ...
 {' ',' ',' ',' ',' ',' ',' ',' '}, ...
 'NumColumns', 2, ...
 'Location', 'best');

lgd.Title.String = 'Experiment          Simulation';
lgd.FontSize = 18;
lgd.FontName = 'Times New Roman';
legend boxon

% --- Annotation ---
annotation('textbox', [0.68 0.12 0.22 0.14], ...
    'String', {'\bfDD6, R_{\epsilon}=-1, (0/0)', '■ 760°C'}, ...
    'Interpreter','tex', ...
    'FontSize',22, ...
    'FontName','Times New Roman', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth',1.2);
hold off;


%% ========================================================================
%% 4. PLOT 2: STRAIN RANGE vs LIFE (NO FIT, TIMES NEW ROMAN, LABEL (a))
%% ========================================================================
strain_range = 2 * [cases.epsilon_a];
X_sim = pred_lives;
X_exp = exp_lives;
Y_vals = strain_range; 

figure('Color','w'); hold on;

% --- PLOT POINTS ONLY ---
hExp = plot(X_exp, Y_vals, 'o', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');

hSim = plot(X_sim, Y_vals, 's', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');

% --- FORMATTING (No Bold, Times New Roman) ---
set(gca, 'XScale', 'log', 'YScale', 'linear', ...
         'FontSize', 22, ...
         'LineWidth', 1.5, ...
         'FontName', 'Times New Roman', ...
         'FontWeight', 'normal'); 

ylim([0.008 0.026]); 
xlim([10 100000]); 
yticks([0.008 0.012 0.016 0.020 0.024]);

xlabel('Cycles to Failure, N_f', ...
    'FontSize', 26, 'FontName', 'Times New Roman', 'FontWeight', 'normal');
ylabel('Strain Range, \Delta\epsilon', ...
    'FontSize', 26, 'FontName', 'Times New Roman', 'FontWeight', 'normal');

% Title removed as requested

grid on; box on;

% --- ADD LABEL (a) Top-Left ---
% Units 'normalized' places it relative to axes (0,0 is bottom-left, 1,1 is top-right)
text(0.02, 0.98, '(a)', ...
    'Units', 'normalized', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');

% --- LEGEND ---
lgd = legend([hExp, hSim], ...
    {'Experiment', 'Simulation'}, ...
    'Location', 'northeast');
lgd.FontSize = 24;
lgd.FontName = 'Times New Roman';
lgd.FontWeight = 'normal';
legend boxon

% --- ANNOTATION (Bottom Left) ---
annotation('textbox', [0.15 0.15 0.22 0.14], ...
    'String', {'DD6, R_{\epsilon}=-1, (0/0)', '■ 760°C'}, ...
    'Interpreter','tex', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth', 1.2);
hold off;

%% ========================================================================
%% 5. PLOT 3: STRAIN-LIFE WITH REGRESSION (Uniform Markers)
%% ========================================================================
figure('Color', 'w'); hold on;

% --- 1. PREPARE DATA ---
X_exp = exp_lives;
Y_exp = 2 * [cases.epsilon_a]; 
X_sim = pred_lives; 
Y_sim = Y_exp; 

% --- 2. PERFORM FITTING (Power Law: Y = A * X^b) ---
log_X = log10(X_exp);
log_Y = log10(Y_exp);
p = polyfit(log_X, log_Y, 1); 
b = p(1); 
A = 10^p(2); 

x_fit_line = logspace(log10(min(X_exp)*0.8), log10(max(X_exp)*1.5), 100);
y_fit_line = A * (x_fit_line).^b;

% --- 3. CALCULATE R-SQUARED (R2) ---
y_fit_at_data = polyval(p, log_X);
y_resid = log_Y - y_fit_at_data;
SS_resid = sum(y_resid.^2);
SS_total = (length(log_Y) - 1) * var(log_Y);
R2 = 1 - (SS_resid / SS_total);

% --- 4. PLOTTING ---
hFit = plot(x_fit_line, y_fit_line, 'k--', 'LineWidth', 2);

hSim = plot(X_sim, Y_sim, 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);

hExp = plot(X_exp, Y_exp, 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

% --- 5. FORMATTING ---
set(gca, 'XScale', 'log', 'YScale', 'log'); 
set(gca, 'LineWidth', 1.5, 'FontSize', 22, 'FontName', 'Times New Roman');

xlabel('Cycles to Failure, N_f', ...
    'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');
ylabel('Strain Range, \Delta\epsilon', ...
    'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');
title(['Coffin-Manson Fit (R^2 = ' sprintf('%.4f', R2) ')'], ...
    'FontWeight', 'bold', 'FontSize', 22, 'FontName', 'Times New Roman');

grid on; box on;
xlim([10 20000]); 
ylim([0.008 0.03]); % Adjusted slightly for 4 cases

% --- 6. LEGEND & ANNOTATION ---
lgd = legend([hExp, hSim, hFit], ...
    {'Experiment', 'Simulation', sprintf('Fit: \\Delta\\epsilon = %.3f N_f^{%.3f}', A, b)}, ...
    'Location', 'southeast');
lgd.FontSize = 18;
lgd.FontName = 'Times New Roman';
legend boxon;

annotation('textbox', [0.65 0.75 0.22 0.1], ...
    'String', {'\bf DD6 Superalloy', '760^{\circ}C'}, ...
    'Interpreter','tex', ...
    'FontSize',18, ...
    'FontName','Times New Roman', ...
    'BackgroundColor','white', ...
    'EdgeColor','black');

hold off;
%% ========================================================================
%% 4. PLOT 2: STRAIN RANGE vs LIFE (WITH POLYNOMIAL FIT)
%% ========================================================================
strain_range = 2 * [cases.epsilon_a];
X_sim = pred_lives;
X_exp = exp_lives;
Y_vals = strain_range; 

figure('Color','w','Position', [900 100 800 800]); hold on;

% --- 1. POLYNOMIAL FIT FOR SIMULATION DATA ---
% We fit log10(Life) vs Strain Range to get a smooth curve on the log scale
valid_idx = ~isnan(X_sim); % Ensure no NaN values interfere
logX = log10(X_sim(valid_idx));
Y_fit_data = Y_vals(valid_idx);

% Fit a 2nd degree polynomial: log10(Nf) = p1*eps^2 + p2*eps + p3
% Or more commonly for S-N: eps = p1*log10(Nf)^2 + ...
p = polyfit(logX, Y_fit_data, 2); 

% Generate points for a smooth dashed line
x_smooth = logspace(log10(min(X_sim)), log10(max(X_sim)), 100);
y_smooth = polyval(p, log10(x_smooth));

% Plot the Blue Dashed Fit Line
hFit = plot(x_smooth, y_smooth, 'b--', 'LineWidth', 2);

% --- 2. PLOT DATA POINTS ---
hExp = plot(X_exp, Y_vals, 'o', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
hSim = plot(X_sim, Y_vals, 's', 'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');

% --- FORMATTING ---
set(gca, 'XScale', 'log', 'YScale', 'linear', ...
         'FontSize', 22, ...
         'LineWidth', 1.5, ...
         'FontName', 'Times New Roman', ...
         'FontWeight', 'normal', 'TickDir', 'in'); 

ylim([0.008 0.026]); 
xlim([10 100000]); 
yticks([0.008 0.012 0.016 0.020 0.024]);
xlabel('Cycles to Failure, N_f', 'FontSize', 26);
ylabel('Strain Range, \Delta\epsilon', 'FontSize', 26);

grid on; box on; axis square;

% Label (c)
text(0.02, 0.98, '(c)', 'Units', 'normalized', 'FontSize', 24, ...
    'FontName', 'Times New Roman', 'VerticalAlignment', 'top');

% --- LEGEND (Updated to include fit line) ---
lgd = legend([hExp, hSim, hFit], ...
    {'Experiment', 'Simulation', 'Sim. Trendline'}, ...
    'Location', 'northeast');
lgd.FontSize = 24;
legend boxon

% --- ANNOTATION (Bottom Left) ---
annotation('textbox', [0.15 0.17 0.3 0.14], ...
    'String', {'DD6, R_{\epsilon}=-1, (0/0)', '■ 760°C'}, ...
    'Interpreter','tex', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'LineWidth', 1.2);

hold off;




%% ========================================================================
%% ========================================================================
%% ========================================================================
%% PLOT: HYSTERESIS LOOPS (ALL CYCLES) - FINAL FIGURE STYLE
%% ========================================================================
figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 850 800]); 
hold on;

% Updated file list and properties to match the reference image style
files = {'cyclic_1_6.csv', 'cyclic_1_8.csv', 'cyclic_2_0.csv', 'cyclic_2_4.csv'};
colors = {'#006400', '#FF0000', '#0000FF', '#00FF00'}; % Green, Red, Blue, Lime
markers = {'none', '>', 'o', 'd'}; 
labels = {'\Delta\epsilon/2 = 0.4%', '\Delta\epsilon/2 = 0.5%', ...
          '\Delta\epsilon/2 = 0.7%', '\Delta\epsilon/2 = 0.9%'};

for k = 1:length(files)
    current_file = files{k};
    
    if isfile(current_file)
        data = readmatrix(current_file);
        strain_all = data(:, 33); 
        stress_all = data(:, 36);
        
        % Plotting with markers and black edges as seen in the reference
        plot(strain_all, stress_all, '-', ...
            'Color', colors{k}, ...
            'LineWidth', 1.5, ...
            'Marker', markers{k}, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', 'none', ... % Hollow markers match your image
            'MarkerEdgeColor', 'k', ...
            'MarkerIndices', round(linspace(1, length(strain_all), 40)), ... 
            'DisplayName', labels{k});
    else
        fprintf('Warning: %s not found.\n', current_file);
    end
end
% 'k:' creates a black dotted line
xline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Vertical dotted line
yline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Horizontal dotted line
% --- FIGURE 4 FORMATTING ---
set(gca, 'LineWidth', 2, ...
         'FontSize', 24, ...
         'FontName', 'Times New Roman', ...
         'TickDir', 'in', ...
         'TickLength', [0.02 0.02], ...
         'Box', 'on', ...
         'XMinorTick', 'on', ...
         'YMinorTick', 'on');

xlabel('Strain, \epsilon', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress, \sigma (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');

axis square;
grid off;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.4);

% Legend - Modified for the requested style
lgd = legend('show', 'Location', 'southeast');
set(lgd, 'FontSize', 22, 'FontName', 'Times New Roman', 'Box', 'on', 'Interpreter', 'tex');




hold off;