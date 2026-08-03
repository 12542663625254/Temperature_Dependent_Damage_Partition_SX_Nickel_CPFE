clc; clear; close all;

%% 1. SETUP: DEFINE MATERIAL CONSTANTS
% =========================================================================
Sg = 0.38; Sc = 0.9*Sg; D0 = 0; Dfc = 0.9; n1 = 0.6;        
Q = 6.97e-19; T = 1033; R_univ = 8.314;       

%% 2. DEFINE FILENAMES
% =========================================================================
% Group A: Pure Fatigue (Dwell = 0/0)
pure_files = {'0_0_R__1.csv', 'mean_0_002.csv', 'mean_0_01.csv', 'mean_0_015.csv'};

% Group B: Dwell Cycles (15/0s or similar)
dwell_files = {'R_1_30_0.csv', 'R_0_67_30_0.csv', 'R_0_30_0.csv', 'R_0_02_30_0.csv'};

% Initialize result arrays
Df_pure = zeros(1, 4);
Df_from_dwell = zeros(1, 4);
D_total_dwell = zeros(1, 4);

%% 3. PROCESSING LOOP
% =========================================================================
for k = 1:4
    % --- CASE 1: Pure Fatigue baseline ---
    % 'false' tells the function to ignore creep logic
    Df_pure(k) = calculate_damage(pure_files{k}, Sg, Sc, D0, Dfc, Q, T, R_univ, false);
    
    % --- CASE 2 & 3: Fatigue vs. Total from Dwell cycle ---
    % We call the function with 'true' to get both components
    [df_dwell, d_total] = calculate_damage_v2(dwell_files{k}, Sg, Sc, D0, Dfc, Q, T, R_univ);
    Df_from_dwell(k) = df_dwell;
    D_total_dwell(k) = d_total;
end

%% 4. VISUALIZATION
% =========================================================================
r_ratios = [-1, -0.67, 0, 0.2]'; 
[r_sort, idx] = sort(r_ratios);

% Sort and Log-protect all three series
plot_pure  = max(Df_pure(idx), 1e-12);
plot_df_dw = max(Df_from_dwell(idx), 1e-12);
plot_total = max(D_total_dwell(idx), 1e-12);

figure('Color','w','Position', [100 100 1100 850]);
hold on;

% Line 1: Pure Fatigue (Blue Circles)
h1 = plot(r_sort, plot_pure, '-ob', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Pure Fatigue ($0/0$ s)');

% Line 2: Fatigue component in Dwell (Green Diamonds)
h2 = plot(r_sort, plot_df_dw, '--dg', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Fatigue Component ($d_f$)');

% Line 3: Total Creep-Fatigue (Red Squares)
h3 = plot(r_sort, plot_total, '-sr', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Total Damage ($d_f + d_c$)');

% Formatting
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylabel('Damage per cycle, $d$', 'Interpreter', 'latex', 'FontSize', 28);
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 28);
ylim([1e-4, 1]); xlim([-1.2, 0.3]);
xticks([-1, -0.67, 0, 0.2]); grid on; box on; axis square;

% Legend & Professional Annotation
legend([h1, h2, h3], 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 18);
annotation('textbox', [0.6 0.15 0.25 0.15], 'String', ...
    {'\textbf{Condition:}', 'Temp: 1033 K', '$\Delta\epsilon = 2.0\%$'}, ...
    'Interpreter', 'latex', 'FontSize', 20, 'BackgroundColor', 'w', 'FitBoxToText','on');

hold off;

%% 5. UPDATED HELPER FUNCTIONS
% =========================================================================

% Standard calculation for single values (Pure Fatigue)
function Df = calculate_damage(file, Sg, Sc, D0, Dfc, Q, T, R_univ, ~)
    [Df, ~] = calculate_damage_v2(file, Sg, Sc, D0, Dfc, Q, T, R_univ);
end

% Detailed calculation for multiple components
function [Df, D_total] = calculate_damage_v2(file, Sg, Sc, D0, Dfc, Q, T, R_univ)
    if ~isfile(file), Df = 1e-12; D_total = 1e-12; return; end
    data = readmatrix(file);
    time = data(:, 1); s_dot_col = data(:, 25);
    eps_a = 0.01; 
    phi = (2.0428e3 + (-2.1056e4 * eps_a) + (-1.2398e7 * eps_a^2)) * exp(-Q / (R_univ * T));
    
    Delta_Sf = 0; Dc = 0;
    for i = 2:length(time)
        dt = time(i) - time(i-1);
        s_dot = s_dot_col(i);
        
        % Logic: Assume Dwell starts after initial loading (t > 5s and s_dot is very small)
        % Adjust '5' or 's_dot' threshold based on your specific cycle data
        if time(i) > 5 && abs(s_dot) < 1e-4 
            Dc = Dc + (1/phi)*(max(s_dot,0))^(0.4)*dt;
        else
            Delta_Sf = Delta_Sf + (s_dot * dt);
        end
    end
    
    % Fatigue calculation
    if Delta_Sf < Sc
        Df = D0 + ((Dfc - D0) / log(1 - Sc/Sg)) * log(1 - Delta_Sf/Sg);
    else
        Df = Dfc;
    end
    
    D_total = Df + Dc;
end

%% 6. VISUALIZATION: LIFE VS STRAIN RATIO
% =========================================================================
% Life is the inverse of damage: N_f = 1 / d
Nf_pure  = 1 ./ plot_pure;
Nf_df_dw = 1 ./ plot_df_dw;
Nf_total = 1 ./ plot_total;

figure('Color','w','Position', [150 150 1100 850]);
hold on; grid on; box on;

% Line 1: Pure Fatigue Life
p1 = plot(r_sort, Nf_pure, '-ob', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Pure Fatigue Life ($0/0$ s)');

% Line 2: Life if only Fatigue occurred during Dwell
p2 = plot(r_sort, Nf_df_dw, '--dg', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Fatigue-only Life (from Dwell)');

% Line 3: Total Combined Life (Creep + Fatigue)
p3 = plot(r_sort, Nf_total, '-sr', 'LineWidth', 3, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Total Creep-Fatigue Life');

% Formatting the Life Axis
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 2);
ylabel('Cycles to Failure, $N_f$', 'Interpreter', 'latex', 'FontSize', 28);
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 28);

% Set limits to capture high fatigue life and low creep-fatigue life
ylim([1, 1e5]); 
xlim([-1.2, 0.3]);
xticks([-1, -0.67, 0, 0.2]);
axis square;

% Legend & Context
legend([p1, p2, p3], 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 18);

title('Life Prediction vs. Strain Ratio', 'FontSize', 24, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

% Add an arrow or annotation to highlight the "Life Reduction"
text(-0.5, 100, '\(\longleftarrow\) Life reduction due to Creep', ...
    'FontSize', 20, 'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'latex');

hold off;