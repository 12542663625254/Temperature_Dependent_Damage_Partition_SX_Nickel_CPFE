clc; clear; close all;

%% 1. SETUP: DEFINE THE FOUR CASES
% =========================================================================
cases(1).filename = '0_0_R__1.csv';   cases(1).epsilon_a = 0.01;  cases(1).exp_life = 30;           
cases(2).filename = 'mean_0_002.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 290;         
cases(3).filename = 'mean_0_01.csv';  cases(3).epsilon_a = 0.01;  cases(3).exp_life = 850;          
cases(4).filename = 'mean_0_015.csv'; cases(4).epsilon_a = 0.01;  cases(4).exp_life = 3557;         

% GLOBAL MATERIAL CONSTANTS
Sg = 0.38; Sc = 0.9*Sg; D0 = 0; Dfc = 0.9; n1 = 0.6;        
Q = 6.97e-19; T = 1033; R_univ = 8.314;       

pred_lives = zeros(1, 4);
Df_vals = zeros(1, 4); 
Dc_vals = zeros(1, 4); 

%% 2. MAIN PROCESSING LOOP
for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    B1 = 2.0428e3 + (-2.1056e4 * eps_a) + (-1.2398e7 * eps_a^2);
    phi = B1 * exp(-Q / (R_univ * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:, 1); entropy_rate_col = data(:, 25); 
        Delta_Sf = 0; Dc = 0;       
        
        for i = 2:length(time_col)
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            % Pure Fatigue Logic for (0/0) cases
            is_fatigue = (t_current > 5 && t_current <= 50); 
            is_creep = false; 
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1 - n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + (s_dot * dt);
            end
        end
        
        % Damage Calculation
        if Delta_Sf < Sc
            Df = D0 + ((Dfc - D0) / log(1 - Sc/Sg)) * log(1 - Delta_Sf/Sg);
        else
            Df = Dfc;
        end
        
        Df_vals(k) = Df;
        Dc_vals(k) = Dc;
        pred_lives(k) = 1 / (Dc + Df);
    end
end

%% 3. VISUALIZATION: REFINED FATIGUE INTERACTION LINE
% =========================================================================
r_ratios = [-1, -0.67, 0, 0.2]'; % R-values
figure('Color','w','Position', [100 100 1000 800], 'Name', 'Fatigue Damage Interaction');
hold on;

% Sort data based on R-ratio
[r_sort, idx] = sort(r_ratios);
Df_sort = Df_vals(idx);

% Prevent Log(0) errors for plotting
Df_plot = max(Df_sort, 1e-12);

% Define custom tick points to match Section 4
y_ticks = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];

% --- PLOT FATIGUE DAMAGE (d_f) ---
h1 = plot(r_sort, Df_plot, '-ob', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Fatigue Damage ($d_f$)');

% --- AXIS FORMATTING ---
set(gca, 'YScale', 'log', 'YColor', 'k');
ylabel('Fatigue damage per cycle, $d_f$', 'Interpreter', 'latex', 'FontSize', 28);
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 28, 'FontName', 'Times New Roman');

% Set limits and ticks
ylim([1e-4, 1]); % Matches Section 4 Y-limit
xlim([-1.5 0.5]);
xticks([-1, -0.67, 0, 0.2]);
yticks(y_ticks);

% --- GLOBAL STYLE ---
grid on; % Enabled to match Section 4
box on; 
axis square;
set(gca, 'FontSize', 22, 'LineWidth', 2, 'FontName', 'Times New Roman', 'TickDir', 'in');
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off');

% Professional Legend
lgd = legend(h1, {'$d_f$ (Fatigue)'}, 'Location', 'northwest', 'Box', 'on', ...
    'EdgeColor', 'k', 'Interpreter', 'latex', 'FontSize', 24);
set(lgd, 'LineWidth', 1.5);

% --- ANNOTATION BOX (Matching Section 4 Position/Style) ---
anno_str = { 'Pure Fatigue Case', ...
             'Dwell = 0/0 s', ...
             'Strain range = 2.0%' };
         
annotation('textbox', [0.62 0.2 0.25 0.15], ...
    'String', anno_str, ...
    'Interpreter', 'latex', ...
    'FontSize', 22, ...
    'FontName', 'Times New Roman', ...
    'BackgroundColor', 'white', ...
    'LineWidth', 1.5, ...
    'EdgeColor', 'k', ...
    'FitBoxToText', 'on');

% Final legend sizing adjustments
currentPos = lgd.Position; 
lgd.Position = [currentPos(1), currentPos(2), currentPos(3)*1.1, currentPos(4)*1.1];

hold off;