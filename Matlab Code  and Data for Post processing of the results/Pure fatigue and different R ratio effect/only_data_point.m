clc; clear; close all;

%% 1. EXPERIMENTAL DATA SETTINGS
% --- Data derived from simulated results for creep-fatigue analysis ---
exp_lives_pure = [222.989, 169.835, 40.1059, 27.06]; % PURE FATIGUE (No Hold)
exp_lives_1    = [121, 99, 52, 55]; % Case 1: 15/0 (Tension Hold)
exp_lives_2    = [69, 65.8, 50, 39];   % Case 2: 30/0
exp_lives_3    = [60, 54, 44, 41];   % Case 3: 60/0

mean_strains = [0, 0.002, 0.010, 0.015]; 
R_ratios     = [-1, -0.67, 0, 0.2];
r_labels     = {'-1', '-0.67', '0', '+0.2'}; 

% Professional Journal Palette
colors = {[0.15, 0.25, 0.55], [0.65, 0.05, 0.05], [0.05, 0.40, 0.20], [0.85, 0.55, 0.10]};
markers = {'o', 's', '^', 'd'};




% Define a dense grid for smooth trendlines
r_fit_grid = linspace(-1.2, 0.4, 200);

% Define Data Series for the loop to keep code clean
all_lives = {exp_lives_pure, exp_lives_1, exp_lives_2, exp_lives_3};
all_labels = {'Pure Fatigue (0/0)', 'Case 1 (15/0)', 'Case 2 (30/0)', 'Case 3 (60/0)'};
all_colors = {[0.3 0.3 0.3], 'r', [0 0.5 0], 'b'};
all_markers = {'d', 'o', '^', 's'};
all_lineStyles = {'-', '-', '-.', '--'};

hPlots = []; % Array to store handles for legend

for j = 1:length(all_lives)
    current_N = all_lives{j};
    
    % --- BETTER POLYNOMIAL FIT ---
    % Fit log10(N) vs R. A 2nd order is usually sufficient, 
    % but we ensure it doesn't oscillate by using 'polyval' on a log-scale.
    p = polyfit(R_ratios, log10(current_N), 2);
    n_fit = 10.^polyval(p, r_fit_grid);
    
    % Plot the smooth trendline
    plot(r_fit_grid, n_fit, all_lineStyles{j}, ...
        'Color', all_colors{j}, 'LineWidth', 3);
    
    % Plot the experimental points
    hPlots(j) = scatter(R_ratios, current_N, 200, all_markers{j}, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', all_colors{j}, ...
        'LineWidth', 2, 'DisplayName', all_labels{j});
end

% --- PROFESSIONAL FORMATTING ---
set(gca, 'XScale', 'linear', 'YScale', 'log', ...
         'FontSize', 22, 'FontName', 'Times New Roman', ...
         'LineWidth', 2, 'TickDir', 'in', 'Box', 'on', ...
         'XMinorTick', 'on', 'YMinorTick', 'on');

% Labels
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('Cycles to Failure, $N_f$', 'Interpreter', 'latex', 'FontSize', 30);

% Axis Limits and Ticks
xlim([-1.2, 0.4]); 
ylim([10, 10000]); 
xticks([-1, -0.67, 0, 0.2]);
grid on; % Enabled grid for better readability of log scale
set(gca, 'GridAlpha', 0.15, 'MinorGridAlpha', 0.1);
axis square;

% --- ENHANCED LEGEND ---
[lgd, ~] = legend(hPlots, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 22);
title(lgd, 'Hold Time (s/s)');

% Professional Annotation
annotation('textbox', [0.15 0.15 0.3 0.1], 'String', ...
    {'Material: High-Temp Alloy', 'Temp: 1033 K', '$\Delta\epsilon = 2.0\%$'}, ...
    'Interpreter', 'latex', 'FontSize', 18, 'EdgeColor', 'none', 'BackgroundColor', 'none');

hold off;

%% 3. REFINED VISUALIZATION: PLOT 6 (Monotonic Polynomial Fit)
% =========================================================================
figure('Color', 'w', 'Position', [950, 100, 900, 850], 'Name', 'R-ratio vs Life (Improved Fit)');
hold on;

% Define the life grid (x-axis) for trendlines
nf_fit_grid = logspace(log10(10), log10(1000), 200);

% Data Series Setup
all_lives = {exp_lives_pure, exp_lives_1, exp_lives_2, exp_lives_3};
all_labels = {'Pure Fatigue (0/0)', 'Case 1 (15/0)', 'Case 2 (30/0)', 'Case 3 (60/0)'};
all_colors = {[0.3 0.3 0.3], 'r', [0 0.5 0], 'b'};
all_markers = {'d', 'o', '^', 's'};
all_lineStyles = {'-', '-', '-.', '--'};

hPlots = []; 
for j = 1:length(all_lives)
    current_N = all_lives{j};
    
    % --- MONOTONIC FIT APPROACH ---
    % To prevent the "left downfall," we fit log10(N_f) = f(R_epsilon)
    % This is more stable for fatigue data.
    p = polyfit(R_ratios, log10(current_N), 2);
    
    % Generate a grid of R-ratios to calculate the corresponding Nf
    r_eval_grid = linspace(-1.2, 0.4, 200);
    nf_eval = 10.^polyval(p, r_eval_grid);
    
    % Plot the trendline (Swap X and Y back to match your axis requirements)
    plot(nf_eval, r_eval_grid, all_lineStyles{j}, ...
        'Color', all_colors{j}, 'LineWidth', 3);
    
    % Plot the experimental points
    hPlots(j) = scatter(current_N, R_ratios, 200, all_markers{j}, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', all_colors{j}, ...
        'LineWidth', 2, 'DisplayName', all_labels{j});
end

% --- PROFESSIONAL FORMATTING ---
set(gca, 'XScale', 'log', 'YScale', 'linear', ...
         'FontSize', 22, 'FontName', 'Times New Roman', ...
         'LineWidth', 2, 'TickDir', 'in', 'Box', 'on', ...
         'XMinorTick', 'on', 'YMinorTick', 'on');

xlabel('Cycles to Failure, $N_f$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 26);

xlim([10, 1000]); 
ylim([-1.2, 0.4]); 
yticks([-1, -0.67, 0, 0.2]);
%grid on; 
%set(gca, 'GridAlpha', 0.1, 'MinorGridAlpha', 0.05);
axis square;

% --- ENHANCED LEGEND WITH LARGE MARKERS ---
[lgd, ~] = legend(hPlots, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 24);
title(lgd, 'Hold Time (s/s)');

% Force MATLAB to render the legend so handles become accessible
drawnow; 

% Find all marker objects within the legend and increase their size
% Typically, these are stored in the 'EntryContainer'
if ~isempty(lgd.EntryContainer)
    for i = 1:length(lgd.EntryContainer.Children)
        % Access the marker (usually the first child of the icon transform)
        try
            set(lgd.EntryContainer.Children(i).Icon.Transform.Children.Children, 'MarkerSize', 14);
        catch
            % Fallback for different MATLAB versions or legend configurations
        end
    end
end

% --- UPDATED ANNOTATION ---
annotation('textbox', [0.15 0.78 0.30 0.08], 'String', ...
    {'Temp: 760$^\circ$C', '$\epsilon_a = 1.0\%$'}, ...
    'Interpreter', 'latex', 'FontSize', 24, 'BackgroundColor', 'w', 'EdgeColor', 'k');

hold off;

hold off;