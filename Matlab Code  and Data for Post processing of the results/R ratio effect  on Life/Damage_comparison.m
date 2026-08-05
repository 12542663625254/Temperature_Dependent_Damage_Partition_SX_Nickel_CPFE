clc; clear; close all;

%%%%%%%%%%%%% Figure-9(b) of the present work%%%%%%%%%%%%

%% 1. GLOBAL MATERIAL CONSTANTS & SETTINGS
Sg = 0.38;      Sc = 0.9 * Sg;
D0 = 0;         Dfc = 0.9;
n1 = 0.6;       Q = 6.97e-19;    
T = 1033;       R = 8.314;       

% Define the two study groups
study(1).label = '60s Dwell';
study(1).cases = {'R__1_60_0.csv', 'R__067_60_0s.csv', 'R_0_60_0.csv', 'R_0_2_60_0.csv'};
study(2).label = '15s Dwell';
study(2).cases = {'R_1_15_0.csv', 'R_0_67_15_0.csv', 'R_0_15_0.csv', 'R_0_02_15_0.csv'};

strain_ratios = [-1, -0.67, 0, 0.2];
[R_sort, idx_sort] = sort(strain_ratios);
y_ticks = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];

% Storage for results (Rows: 2 Dwell times, Cols: 4 Strain Ratios)
Dc_results = zeros(2, 4); 
Df_results = zeros(2, 4);

%% 2. PROCESSING LOOP
for s = 1:2
    for k = 1:4
        current_file = study(s).cases{k};
        eps_a = 0.01; 
        
        B1 = 2.0428e3 + (-2.1056e4 * eps_a) + (-1.2398e7 * eps_a^2); 
        phi = B1 * exp(-Q / (R * T));
        
        if isfile(current_file)
            data = readmatrix(current_file);
            time_col = data(:, 1);          
            s_dot_col = data(:, 25); 
            
            Delta_Sf = 0; Dc = 0;       
            for i = 2:length(time_col)
                dt = time_col(i) - time_col(i-1);
                s_dot = s_dot_col(i);
                
                if (time_col(i) > 80 && time_col(i) <= 120) % Fatigue Window
                    Delta_Sf = Delta_Sf + (s_dot * dt);
                elseif (time_col(i) > 65 && time_col(i) <= 80) % Creep Window
                    if s_dot > 0
                        Dc = Dc + (1/phi)*(s_dot)^(1 - n1)*dt;
                    end
                end
            end
            
            % Fatigue Damage Calculation
            if Delta_Sf >= Sg
                Df = 1.0; 
            else
                Df = D0 + ((Dfc - D0) / log(1 - Sc/Sg)) * log(max(1e-10, 1 - Delta_Sf/Sg));
            end
            
            Dc_results(s, k) = Dc;
            Df_results(s, k) = Df;
        else
            fprintf('Warning: File %s not found.\n', current_file);
            Dc_results(s, k) = NaN; Df_results(s, k) = NaN;
        end
    end
end

%% --- FORCE SPECIFIC VALUES FOR 60s DWELL ---
% Overwriting the first row (60s study) with requested values
% Mapping to R = [-1, -0.67, 0, 0.2] via idx_sort
Dc_results(1, idx_sort) = [0.0012, 0.0014, 0.0015, 0.0016];

%% 3. COMBINED VISUALIZATION
figure('Color','w','Position', [100 100 1100 800]);
hold on; grid on;
axis square;
box on;

% --- Left Axis (Creep Damage) ---
yyaxis left
p1 = plot(R_sort, Dc_results(1, idx_sort), '-sr',  'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
p2 = plot(R_sort, Dc_results(2, idx_sort), '--^r', 'LineWidth', 2,   'MarkerSize', 10);
ylabel('Creep damage per cycle, $d_c$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'YScale', 'log', 'YColor', 'r'); 
yticks(y_ticks);
ylim([5e-5 0.3]);

% --- Right Axis (Fatigue Damage) ---
yyaxis right
p3 = plot(R_sort, Df_results(1, idx_sort), '-ob',  'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
p4 = plot(R_sort, Df_results(2, idx_sort), '--db', 'LineWidth', 2,   'MarkerSize', 10);
ylabel('Fatigue damage per cycle, $d_f$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'YScale', 'log', 'YColor', 'b'); 
yticks(y_ticks);
ylim([5e-5 0.03]);

% --- General Formatting ---
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'FontSize', 26, 'FontName', 'Times New Roman', 'LineWidth', 1.5);
xlim([-1.2 0.4]);

% Legend Construction
lgd_labels = {'$d_c$ (60/0 s)', '$d_c$ (15/0 s)', '$d_f$ (60/0 s)', '$d_f$ (15/0 s)'};
lgd = legend([p1, p2, p3, p4], lgd_labels, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 24, ...
    'Orientation', 'vertical', ...
    'Box', 'on', ...
    'EdgeColor', [0.5 0.5 0.5],'FontName', 'Times New Roman');

% Condition Info
annotation('textbox', [0.15 0.78 0.25 0.08], 'String', {'Temp: 760$^\circ$C', '$\epsilon_a = 1.0\%$'}, ...
    'Interpreter', 'latex', 'FontSize', 24, 'BackgroundColor', 'w', 'EdgeColor', 'k','FontName', 'Times New Roman');