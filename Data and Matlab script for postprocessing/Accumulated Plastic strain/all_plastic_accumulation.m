clc; clear; close all;

%%%%%%% Figure-8 of the present work%%%%%%%%%

%% 1. SETUP & DATA LOADING
% File names for DD6 Superalloy
files = {'0_0_2_4.csv', '30_0_2_4.csv', '60_0_2_4.csv', 'CFI_2_4_120_0.csv', '240_0_2_4.csv', '300_0_2_4_new.csv'};
labels = {'0/0 Fatigue', '30/0 Hold', '60/0 Hold', '120/0 Hold', '240/0 Hold', '300/0 Hold'};

% Experimental Cycle Times (Verify these durations)
T_cycs = [48, 78, 108, 168, 288, 348]; 
t_start = 12;

% Professional Colors and Markers
colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980; % Red-Orange
    0.4660, 0.6740, 0.1880; % Green
    0.9290, 0.6940, 0.1250; % Yellow-Gold
    0.4940, 0.1840, 0.5560; % Purple
    1, 0, 1;                % Magenta
];
marker_styles = {'o', 's', '^', 'v', 'd', 'h'};

%% 2. DATA PROCESSING LOOP (Creates the 'ds' variable)
ds = struct(); 
for i = 1:length(files)
    if ~isfile(files{i})
        fprintf('Warning: %s not found. Skipping...\n', files{i});
        continue;
    end
    
    raw_data = readmatrix(files{i});
    ds(i).time = raw_data(:, 1);
    % Using Column 5 as per your previous calculation logic
    ds(i).s_dot = raw_data(:, 5); 
    
    % Calculations: Cycles and Cumulative Sum
    ds(i).cycles = (ds(i).time - t_start) / T_cycs(i) + 1;
    ds(i).s_acc = cumtrapz(ds(i).time, ds(i).s_dot);
end

%% 3. FIGURE 2: CUMULATIVE ACCUMULATION (Cycle-domain)
figure('Color', 'w', 'Name', 'Cumulative Comparison', 'Position', [150 150 1000 800]);
hold on;

target_cycles = 1:1:10; 

for i = 1:length(ds)
    if isempty(ds(i).time), continue; end
    
    % Locate indices for markers exactly at integer cycles
    marker_idx = zeros(1, length(target_cycles));
    for j = 1:length(target_cycles)
        [~, idx] = min(abs(ds(i).cycles - target_cycles(j)));
        marker_idx(j) = idx;
    end
    
    % Main Plotting Command
    plot(ds(i).cycles, ds(i).s_acc, ...
         'LineStyle', '-', 'LineWidth', 3, 'Color', colors(i,:), ...
         'Marker', marker_styles{i}, ...
         'MarkerSize', 12, ...
         'MarkerEdgeColor', colors(i,:), ... 
         'MarkerFaceColor', 'w', ...
         'MarkerIndices', marker_idx, ... 
         'DisplayName', labels{i});
end

% --- PLOT FORMATTING ---
xlim([1 10]); 
% Set ylim to [0 2] or use [0 0.002] depending on your actual data units
ylim([0 2]); 
grid off; box on; axis square;

set(gca, 'FontSize', 22, 'LineWidth', 2, 'FontName', 'Times New Roman', 'TickDir', 'in');
set(gca, 'XTick', 1:10); 

xlabel('Number of cycles, $N$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Cumulated inelastic strain (\%)', 'Interpreter', 'latex', 'FontSize', 26);

% Legend configuration
legend('Location', 'southeast', 'FontSize', 18, 'Interpreter', 'none'); 

fprintf('Success: Figure 2 generated with all 6 cases and matching markers.\n');