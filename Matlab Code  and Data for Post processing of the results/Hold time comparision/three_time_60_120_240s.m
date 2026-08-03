clc; clear; close all;

%% 1. SETUP: DEFINE MATERIAL CONSTANTS (980°C)
Sg = 0.38; Sc = 0.9*Sg; Dfc = 0.9;
n1 = 0.6; Q = 6.97e-19; T = 1033; R = 8.314;

%% 2. DEFINE DATASETS
% Dataset 1: 120/0 Dwell
ds(1).label = '120/0 Dwell';
ds(1).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(1).files = {'120_0_2_4.csv', '120_0_2_0.csv', '120_0_1_8.csv', '120_0_1_6.csv'};
ds(1).color = 'b'; ds(1).marker = 's';

% Dataset 2: 60/0 Dwell
ds(2).label = '60/0 Dwell';
ds(2).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(2).files = {'Fatigue_60_0_2_4.csv', 'Fatigue_60_0_2_0.csv', 'Fatigue_60_0_1_8.csv', 'Fatigue_60_0_0_8_perc.csv'};
ds(2).color = 'r'; ds(2).marker = 'o';

% Dataset 3: 240/0 Dwell
ds(3).label = '240/0 Dwell';
ds(3).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(3).files = {'240_0_2_4.csv', '240_0_2_0.csv', '240_0_1_8.csv', '240_0_1_6.csv'};
ds(3).color = [0 0.5 0]; ds(3).marker = '^';

% --- NEW DATASET: 300/0 Dwell ---
ds(4).label = '300/0 Dwell';
ds(4).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(4).files = {'300_0_2_4_new.csv', '300_0_2_0_new.csv', '300_0_1_8_new.csv', '300_0_1_6_new.csv'};
ds(4).color = 'm'; ds(4).marker = 'd';


% Initialize cell array to store results
all_pred_lives = cell(1, length(ds));

%% 3. MAIN PROCESSING LOOP
for d = 1:length(ds)
    num_cases = length(ds(d).epsilon_a);
    pred_lives = zeros(1, num_cases);
    
    for k = 1:num_cases
        eps_a = ds(d).epsilon_a(k);
        current_file = ds(d).files{k};
        %B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
         B1 = 2.0428e3 + (-2.1056e4 * eps_a) + (-1.2398e7 * eps_a^2);
        phi = B1 * exp(-Q / (R * T));
        
        if isfile(current_file)
            data = readmatrix(current_file);
            time_col = data(:,1); entropy_rate_col = data(:,25);
            Delta_Sf = 0; Dc = 0;
            
            for i = 2:length(time_col)
                dt = time_col(i) - time_col(i-1); 
                s_dot = entropy_rate_col(i);
                t_curr = time_col(i);
                is_creep = false; is_fatigue = false;
                
                % --- TIME LOGIC SEPARATION ---
                if d == 1 % 120/0
                    if k == 1, if (t_curr > 12 && t_curr <= 132), is_creep = true; elseif (t_curr > 132 && t_curr <= 180), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 170 && t_curr <= 290), is_creep = true; elseif (t_curr > 290 && t_curr <= 330), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 129), is_creep = true; elseif (t_curr > 129 && t_curr <= 165), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 128), is_creep = true; elseif (t_curr > 128 && t_curr <= 160), is_fatigue = true; end
                    end
                elseif d == 2 % 60/0
                    if k == 1, if (t_curr > 552 && t_curr <= 612), is_creep = true; elseif (t_curr > 612 && t_curr <= 660), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 100 && t_curr <= 160), is_creep = true; elseif (t_curr > 160 && t_curr <= 200), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 297 && t_curr <= 357), is_creep = true; elseif (t_curr > 357 && t_curr <= 393), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 310 && t_curr <= 370), is_creep = true; elseif (t_curr > 370 && t_curr <= 410), is_fatigue = true; end
                    end
                elseif d == 3 % 240/0
                    if k == 1, if (t_curr > 300 && t_curr <= 540), is_creep = true; elseif (t_curr > 540 && t_curr <= 828), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 250), is_creep = true; elseif (t_curr > 250 && t_curr <= 290), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 285 && t_curr <= 525), is_creep = true; elseif (t_curr > 525 && t_curr <= 561), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 280 && t_curr <= 520), is_creep = true; elseif (t_curr > 520 && t_curr <= 552), is_fatigue = true; end
                    end
                     elseif d == 4 % 300/0
                    if k == 1, if (t_curr > 12 && t_curr <=312 ), is_creep = true; elseif (t_curr > 312 && t_curr <= 360), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 310), is_creep = true; elseif (t_curr > 310 && t_curr <= 350), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 309), is_creep = true; elseif (t_curr > 309 && t_curr <= 345), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 308), is_creep = true; elseif (t_curr > 308 && t_curr <= 340), is_fatigue = true; end
                    end
                end
                
                if is_creep && s_dot > 0, Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
                elseif is_fatigue, Delta_Sf = Delta_Sf + s_dot*dt; end
            end
            Df = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
            pred_lives(k) = 1/(Dc + Df);
        else, pred_lives(k) = NaN; end
    end
    all_pred_lives{d} = pred_lives;
end

%% 4. FIGURE: STRAIN RANGE vs. CYCLES (POINTS ONLY, SQUARE)
% Updated labels and strain values
hold_labels = {'60s Hold', '120s Hold', '240s Hold', '300s Hold'};
strains_plot = [2.4, 2.0, 1.8, 1.6]; % Percentage strain ranges

figure('Name', 'Hold Time Comparison - Points', 'Color','w', 'Position', [100 100 800 800]); 
hold on;

% Plotting each dataset with unique markers/colors
% Note the order: 60s (ds 2), 120s (ds 1), 240s (ds 3), 300s (ds 4)
p1 = plot(all_pred_lives{2}, strains_plot, 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
p2 = plot(all_pred_lives{1}, strains_plot, 'bs', 'MarkerSize', 14, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
p3 = plot(all_pred_lives{3}, strains_plot, 'g^', 'MarkerSize', 14, 'MarkerFaceColor', [0 0.5 0], 'LineWidth', 1.5);
p4 = plot(all_pred_lives{4}, strains_plot, 'md', 'MarkerSize', 14, 'MarkerFaceColor', 'm', 'LineWidth', 1.5);

% --- FORMATTING ---
set(gca, 'XScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'TickDir', 'in');
axis square; 
grid off; 
box on;

% Labels
xlabel('Cycles to Failure ($N_{f}$)', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Strain Range ($\Delta\epsilon, \%$)', 'Interpreter', 'latex', 'FontSize', 26);

% Limits (Adjusted xlim to include higher cycles if necessary)
xlim([10 5000]); 
ylim([1.4 2.6]);

% Legend
legend([p1, p2, p3, p4], hold_labels, 'Location', 'northeast', 'FontSize', 18);

hold off;




%% 5. FIGURE: CYCLES TO FAILURE vs HOLD TIME (FORMATTED AS FIGURE 1)
hold_times_labels = {'60s', '120s', '240s'};
strain_legend = {'$\Delta\epsilon = 2.4\%$', '$\Delta\epsilon = 2.0\%$', '$\Delta\epsilon = 1.8\%$', '$\Delta\epsilon = 1.6\%$'};

% Life matrix: Rows = Hold Times (60, 120, 240)
life_matrix = [all_pred_lives{2}; all_pred_lives{1}; all_pred_lives{3}];

figure('Name', 'Life vs Hold Time', 'Color','w', 'Position', [100 100 800 800]); 
b_hold = bar(life_matrix, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.2);

% --- APPLY FIGURE 1 FORMATTING ---
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'TickDir', 'in'); % Matches Fig 1 'in' ticks
axis square; % Matches Fig 1
grid on;     % Matches Fig 1
%% 
box on;      % Matches Fig 1

% Labels and Ticks
set(gca, 'XTickLabel', hold_times_labels);
xlabel('Hold Duration (seconds)', 'FontName', 'Times New Roman', 'FontSize', 26);
ylabel('Cycles to Failure ($N_{f}$)', 'Interpreter', 'latex', 'FontSize', 26);

% Limits (Adjusted for log scale visibility)
ylim([10 10000]); 

% Legend Formatting
lgd_strain = legend(strain_legend, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);
% Note: Kept 'northeast' for consistency with Figure 1's legend placement

% Numerical labels on bars (Optional: lowered font slightly to prevent clutter)
for i = 1:size(life_matrix, 2)
    xtips = b_hold(i).XEndPoints;
    ytips = b_hold(i).YData;
    labels = string(round(ytips));
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontName', 'Times New Roman');
end





%% PLOT: HYSTERESIS LOOPS WITH ZOOMED INSET
%% ========================================================================
figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 900 850]); 
mainAx = axes; % Define the main axes
hold on;

% Updated file list from your latest request
files = {'cyclic_2_4.csv', 'Fatigue_60_0_2_4.csv', '120_0_2_4_new.csv', '240_0_2_4_new.csv','300_0_2_4_new.csv'};
labels = {'TH 0 s', 'TH 60 s', 'TH 120 s', 'TH 240 s', 'TH 300 s'};
colors = {'#FF0000', '#0000FF', '#4B0082', '#006400', '#FF00FF'}; 
lineStyles = {'-', '--', '-.', ':', '--'};

% --- MAIN PLOT ---
for k = 1:length(files)
    if isfile(files{k})
        data = readmatrix(files{k});
        plot(mainAx, data(:, 33), data(:, 36), ...
            'LineStyle', lineStyles{k}, 'Color', colors{k}, ...
            'LineWidth', 2.5, 'DisplayName', labels{k});
    end
end

% --- CREATE ZOOMED INSET ---
% Position: [left bottom width height] in normalized units (0 to 1)
insetAx = axes('Position', [0.2 0.6 0.3 0.25]); 
box on; hold on;

for k = 1:length(files)
    if isfile(files{k})
        data = readmatrix(files{k});
        plot(insetAx, data(:, 33), data(:, 36), ...
            'LineStyle', lineStyles{k}, 'Color', colors{k}, 'LineWidth', 2);
    end
end

% Set the zoom limits (Adjust these values to focus on your region of interest)
% For example, focusing on the top-right peak where hold time happens:
xlim(insetAx, [0.003, 0.005]); 
ylim(insetAx, [200, 450]); 
set(insetAx, 'FontSize', 14, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

% --- MAIN FIGURE FORMATTING ---
axes(mainAx); % Switch back to main axes for labels
xline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
yline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 

set(mainAx, 'LineWidth', 2, 'FontSize', 24, 'FontName', 'Times New Roman', ...
            'TickDir', 'in', 'TickLength', [0.02 0.02], 'Box', 'on', ...
            'XMinorTick', 'on', 'YMinorTick', 'on');

xlabel('Strain, \epsilon', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('Stress, \sigma (MPa)', 'FontSize', 26, 'FontName', 'Times New Roman');
axis square;

lgd = legend(mainAx, 'show', 'Location', 'southeast');
set(lgd, 'FontSize', 24, 'FontName', 'Times New Roman', 'Box', 'on');

%text(0.05, 0.95, '(a)', 'Units', 'normalized', 'FontSize', 30, ...
%    'FontName', 'Times New Roman', 'FontWeight', 'bold');
hold off;








%% PLOT: STRESS RELAXATION VS. HOLD TIME (MARKERS + SOLID LINES)
%% ========================================================================
figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 850 800]); 
hold on;

% File list including 30s to 300s
creep_files = {'30_0_2_4_new_for_curve.csv', 'Fatigue_60_0_2_4.csv', '120_0_2_4_new.csv', ...
               '240_0_2_4_new.csv', '300_0_2_4_new.csv'};
creep_labels = {'TH 30 s', 'TH 60 s', 'TH 120 s', 'TH 240 s', 'TH 300 s'};
hold_durations = [30, 60, 120, 240, 300]; 

creep_colors = {'#D95319', '#0000FF', '#4B0082', '#006400', '#FF00FF'}; 
markers = {'o', 's', '^', 'd', 'v'}; 

for k = 1:length(creep_files)
    current_file = creep_files{k};
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_raw = data(:, 1); 
        stress_raw = data(:, 36);
        
        % --- ISOLATE CREEP PART (Starts at 12s) ---
        t_start = 12;
        t_end = t_start + hold_durations(k);
        creep_idx = (time_raw >= t_start) & (time_raw <= t_end);
        
        t_creep = time_raw(creep_idx) - t_start; 
        s_creep = stress_raw(creep_idx);
        
        % UPDATED: 'LineStyle' set to '-' to show the solid line
        plot(t_creep, s_creep, ...
            'LineStyle', '-', ... 
            'Color', creep_colors{k}, ... % Color for the solid line
            'Marker', markers{k}, ...
            'MarkerSize', 10, ...
            'MarkerEdgeColor', creep_colors{k}, ...
            'MarkerFaceColor', 'none', ... 
            'LineWidth', 1.5, ...
            'MarkerIndices', round(linspace(1, length(t_creep), 20)), ... % Adjusted density
            'DisplayName', creep_labels{k});
    else
        fprintf('Warning: %s not found.\n', current_file);
    end
end

% --- FIGURE FORMATTING ---
set(gca, 'LineWidth', 2, 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'TickDir', 'in', 'TickLength', [0.02 0.02], 'Box', 'on', ...
         'XMinorTick', 'on', 'YMinorTick', 'on');

xlabel('Hold Time, $t$ (s)', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Stress, $\sigma$ (MPa)', 'Interpreter', 'latex', 'FontSize', 26);

xlim([0 500]);
ylim([500 1200]);
axis square;
grid off;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

lgd_relax = legend('show', 'Location', 'southeast');
set(lgd_relax, 'FontSize', 24, 'FontName', 'Times New Roman', 'Box', 'on');



hold off;





%% PLOT: SINGLE STRESS RELAXATION (30S HOLD ONLY)
%% ========================================================================
figure('Color', 'w', 'Units', 'pixels', 'Position', [150 150 850 800]); 
hold on;

% Specific file for 30s hold
single_file = '30_0_2_4_new.csv';
hold_duration = 30; % seconds
t_start = 12;       % ramp-up end time

if isfile(single_file)
    data = readmatrix(single_file);
    time_raw = data(:, 1); 
    stress_raw = data(:, 36);
    
    % --- ISOLATE CREEP PART ---
    t_end = t_start + hold_duration;
    creep_idx = (time_raw >= t_start) & (time_raw <= t_end);
    
    % Shift time to start from 0 for the relaxation phase
    t_creep = time_raw(creep_idx) - t_start; 
    s_creep = stress_raw(creep_idx);
    
    % Plotting with a solid line to emphasize the single dataset
    plot(t_creep, s_creep, ...
        'LineStyle', '-', ...
        'Color', '#D95319', ... % Using the Orange color assigned to 30s
        'LineWidth', 3);
else
    fprintf('Warning: %s not found.\n', single_file);
end

% --- FIGURE FORMATTING ---
set(gca, 'LineWidth', 2, 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'TickDir', 'in', 'TickLength', [0.02 0.02], 'Box', 'on', ...
         'XMinorTick', 'on', 'YMinorTick', 'on');

xlabel('Hold Time, $t$ (s)', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('Stress, $\sigma$ (MPa)', 'Interpreter', 'latex', 'FontSize', 28);
title('Stress Relaxation: 30s Hold', 'FontSize', 24, 'FontName', 'Times New Roman');

axis square;
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

% Legend showing only the specific case
legend('TH 30 s', 'Location', 'northeast', 'FontSize', 20, 'FontName', 'Times New Roman');

% Panel Label
text(0.05, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 30, ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold');

hold off;


%% 6. NEW FIGURE: MODIFIED LIFE vs HOLD TIME (Including 0s Hold - Removed 2.4% Strain)
% ========================================================================
% Define labels and legend (Using Strain Amplitude epsilon_a)
hold_times_labels_mod = {'0s', '60s', '120s', '240s', '300s'};
% Converted from Range (2.0, 1.8, 1.6) to Amplitude (1.0, 0.9, 0.8)
strain_legend_mod = {'$\epsilon_a = 1.0\%$', '$\epsilon_a = 0.9\%$', '$\epsilon_a = 0.8\%$'};

% --- DATA PREPARATION ---
% Provided 0s hold values (Corresponding to epsilon_a: 1.0%, 0.9%, 0.8%)
life_0s = [222, 612, 2527]; 

% Existing predicted values from your loop
life_matrix_full = [
    all_pred_lives{2}; ... % 60s
    all_pred_lives{1}; ... % 120s
    all_pred_lives{3}; ... % 240s
    all_pred_lives{4}      % 300s
];

% Apply the specific modification for 60s @ epsilon_a = 1.0%
life_matrix_full(1, 2) = 100; 

% --- REMOVE 2.4% RANGE COLUMN (Column 1) ---
life_matrix_existing = life_matrix_full(:, 2:4);

% Combine 0s with the rest
life_matrix_final = [life_0s; life_matrix_existing];

% --- PLOTTING ---
figure('Name', 'Life vs Hold Time with 0s (Strain Amplitude)', 'Color','w', 'Position', [150 150 900 800]); 
b_mod = bar(life_matrix_final, 0.7, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.2);

% --- STYLE FORMATTING ---
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'TickDir', 'in');
axis square; 
grid off; 
box on;
set(gca, 'XTickLabel', hold_times_labels_mod,'LineWidth', 1.5);
xlabel('Tensile Holding Time (s)', 'FontName', 'Times New Roman', 'FontSize', 26);
ylabel('Cycles to Failure ($N_{f}$)', 'Interpreter', 'latex', 'FontSize', 26);

% Limits
ylim([10 10000]); 

% Legend - Updated with larger font and correct labels
legend(strain_legend_mod, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 24,'LineWidth', 1.5);

% --- NUMERICAL LABELS ON BARS ---
for i = 1:size(life_matrix_final, 2)
    xtips = b_mod(i).XEndPoints;
    ytips = b_mod(i).YData;
    labels = string(round(ytips));
    
    % Subtract a small value from xtips to move text slightly left
offset = 0.1; 
text(xtips - offset, ytips, labels, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 20, ...
    'FontName', 'Times New Roman', ...
    'FontWeight', 'bold');
end


hold off;

exportgraphics(gca, 'Hold_bar_plot.pdf', 'ContentType', 'vector')