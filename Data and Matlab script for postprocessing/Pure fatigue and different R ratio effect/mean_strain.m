clc; clear; close all;

%% 1. SETTINGS & FILENAMES
files = {'mean_strain_0.csv', 'mean_strain_0_002.csv', ...
         'mean_strain_0_01.csv', 'mean_strain_0_015.csv'};                    
starts = [1900, 1500, 1001, 973]; 
mean_strains = [0, 0.002, 0.010, 0.015]; 
eps_a_vals = [0.010, 0.002, 0.010, 0.015]; 

% Standardized labels for legends
labels = {'\epsilon_{mean} = 0', '\epsilon_{mean} = 0.002', '\epsilon_{mean} = 0.010', '\epsilon_{mean} = 0.015'};

% --- CALCULATE R-RATIOS ---
r_labels = {'-1', '-0.67', '0', '0.2'}; 

%% 2. GLOBAL MATERIAL CONSTANTS
n1 = 0.6; Q = 6.97e-19; T = 1033; R_const = 8.314;       
Sg = 0.055; Sc = 0.01268; Dfc = 0.055;  

%% 3. DATA PROCESSING & LIFE CALCULATION
processed_data = cell(1, 4);
pred_lives = zeros(1, 4);
for k = 1:4
    if isfile(files{k})
        raw = readmatrix(files{k});
        time_col = raw(:,1); s_dot_col = raw(:,25);     
        stress_col = raw(:,36); strain_col = raw(:,33);    
        
        processed_data{k}.t_full = time_col;
        processed_data{k}.e_full = strain_col;
        
        idx_start = find(stress_col >= starts(k), 1, 'first');
        if isempty(idx_start), idx_start = 1; end
        
        if k == 4
            [~, max_idx_rel] = max(stress_col(idx_start:end));
            idx_end = idx_start + max_idx_rel - 1;
        else
            idx_end = length(stress_col);
        end
        
        processed_data{k}.x = strain_col(idx_start:idx_end);
        processed_data{k}.y = stress_col(idx_start:idx_end);
        
        eps_a = eps_a_vals(k);
        B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
        phi = B1 * exp(-Q / (R_const * T));
        
        Delta_Sf = 0; Dc = 0;
        for j = 2:length(time_col)
            t_curr = time_col(j);
            dt = time_col(j) - time_col(j-1);
            s_dot = s_dot_col(j);
            
            if k == 1
                is_creep = (t_curr > 110 && t_curr <= 170);
                is_fatigue = (t_curr > 170 && t_curr <= 210);
            elseif k == 2
                is_creep = (t_curr > 312 && t_curr <= 372);
                is_fatigue = (t_curr > 372 && t_curr <= 412);
            elseif k == 3
                is_creep = (t_curr > 20 && t_curr <= 80);
                is_fatigue = (t_curr > 80 && t_curr <= 120);
            elseif k == 4
                is_creep = (t_curr > 25 && t_curr <= 85);
                is_fatigue = (t_curr > 85 && t_curr <= 125);
            end
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + s_dot*dt;
            end
        end
        Df = (Dfc/log(1-Sc/Sg))*log(1-min(Delta_Sf, Sg-1e-5)/Sg);
        pred_lives(k) = 1/(Dc + Df);
    end
end

%% 4. PROFESSIONAL JOURNAL PALETTE & COMMON SETTINGS
journal_colors = {[0.15, 0.25, 0.55], [0.65, 0.05, 0.05], [0.05, 0.40, 0.20], [0.85, 0.55, 0.10]};
bar_width = 0.5;
mkr = {'o', 's', '^', 'd'}; styles = {'-', '--', ':', '-.'};
markers = {'o', 's', '^', 'd'};
styles  = {'-', '--', ':', '-.'};
%% 3. PLOT STRESS–STRAIN HYSTERESIS
%% ============================================================
figure('Color','w','Position',[100 100 850 850]);
hold on;

for k = 1:4
    
    if isfile(files{k})
        
        raw = readmatrix(files{k});
        
        stress_col = raw(:,36);   % Stress column
        strain_col = raw(:,33);   % Strain column
        
        % Find starting index
        idx_start = find(stress_col >= starts(k), 1, 'first');
        if isempty(idx_start)
            idx_start = 1;
        end
        
        % Use remaining data
        strain_plot = strain_col(idx_start:end);
        stress_plot = stress_col(idx_start:end);
        
        plot(strain_plot, stress_plot, ...
    'Color', journal_colors{k}, ...
    'LineStyle', styles{k}, ...
    'LineWidth', 2.2, ...
    'Marker', markers{k}, ...
    'MarkerSize', 9, ...
    'MarkerFaceColor', 'w', ...
    'MarkerIndices', round(linspace(1, length(strain_plot), 40)));
        
    else
        warning('%s not found.', files{k});
    end
    
end

%% ============================================================
%% 4. FORMATTING
%% ============================================================
set(gca, 'FontSize', 24, ...
         'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, ...
         'TickDir', 'in');

xlabel('Strain, $\epsilon$', ...
       'Interpreter','latex', ...
       'FontSize', 26);

ylabel('Stress, $\sigma$ (MPa)', ...
       'Interpreter','latex', ...
       'FontSize', 26);

legend(labels, ...
       'Location','northeast', ...
       'FontSize',22, ...
      'FontName','Times New Roman');
%yticks([-1200 -600 0 50]);
xlim([-0.03 0.03]);
grid on;
box on;
axis square;
text(0.02, 0.98, '(a)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold off;
%% 6. PLOT (b): LIFE vs MEAN STRAIN (LOG-BAR)
figure('Color','w','Position', [100 100 850 850]); hold on;
b = bar(1:4, pred_lives, bar_width, 'EdgeColor', 'k', 'LineWidth', 1.5, 'BaseValue', 1);
b.FaceColor = 'flat';
for i = 1:4, b.CData(i,:) = journal_colors{i}; end
for i = 1:4
    text(i, pred_lives(i)*1.15, sprintf('%.3f', pred_lives(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 20, 'FontName', 'Times New Roman');
end
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Mean Strain, $\epsilon_{mean}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid on; box on; axis square; 
set(gca, 'XTick', 1:4, 'XTickLabel', string(mean_strains));
xlim([0.5, 4.5]); ylim([1, max(pred_lives)*20]);
text(0.02, 0.98, '(c)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%% 7. PLOT (c): STRAIN vs TIME
figure('Color','w','Position', [100 100 850 850]); hold on;
h_time = [];
for i = 1:4
    if ~isempty(processed_data{i})
        h_time(i) = plot(processed_data{i}.t_full, processed_data{i}.e_full, ...
            'Color', journal_colors{i}, 'LineStyle', styles{i}, 'LineWidth', 2.0);
    end
end
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Time, $t$ (s)', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Strain, $\epsilon$', 'Interpreter', 'latex', 'FontSize', 26);
grid off; box on; axis square; xlim([0 210]); 
text(0.02, 0.98, '(b)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
axis square;
lgd3 = legend(h_time, labels, 'Location', 'northeast');
lgd3.FontSize = 24; lgd3.FontName = 'Times New Roman'; legend boxon;

%% 8. PLOT (d): LIFE vs R-RATclc; clear; close all;
%% 1. SETTINGS & FILENAMES
files = {'Fatigue.csv', '0_002_60_0.csv', ...
         'Mean_strain_010_1_percentage.csv', '0_015_60_0.csv'};                    
starts = [1900, 1500, 1810, 2500]; 
mean_strains = [0, 0.002, 0.010, 0.015]; 
eps_a_vals = [0.010, 0.002, 0.010, 0.015]; 
labels = {'\epsilon_{m} = 0', '\epsilon_{m} = 0.002', '\epsilon_{m} = 0.010', '\epsilon_{m} = 0.015'};
r_vals = [-1, -0.67, 0, 0.2]; % Numerical values for fitting
r_labels = {'-1', '-0.67', '0', '0.2'}; 

%% 2. GLOBAL MATERIAL CONSTANTS
n1 = 0.6; Q = 6.97e-19; T = 1033; R_const = 8.314;       
Sg = 0.055; Sc = 0.01268; Dfc = 0.055;  

%% 3. DATA PROCESSING & LIFE CALCULATION
processed_data = cell(1, 4);
pred_lives = zeros(1, 4);
for k = 1:4
    if isfile(files{k})
        raw = readmatrix(files{k});
        time_col = raw(:,1); s_dot_col = raw(:,25);     
        stress_col = raw(:,36); strain_col = raw(:,33);    
        
        processed_data{k}.t_full = time_col;
        processed_data{k}.e_full = strain_col;
        
        idx_start = find(stress_col >= starts(k), 1, 'first');
        if isempty(idx_start), idx_start = 1; end
        
        if k == 4
            [~, max_idx_rel] = max(stress_col(idx_start:end));
            idx_end = idx_start + max_idx_rel - 1;
        else
            idx_end = length(stress_col);
        end
        
        processed_data{k}.x = strain_col(idx_start:idx_end);
        processed_data{k}.y = stress_col(idx_start:idx_end);
        
        eps_a = eps_a_vals(k);
        B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
        phi = B1 * exp(-Q / (R_const * T));
        
        Delta_Sf = 0; Dc = 0;
        for j = 2:length(time_col)
            t_curr = time_col(j);
            dt = time_col(j) - time_col(j-1);
            s_dot = s_dot_col(j);
            
            if k == 1
                is_creep = (t_curr > 110 && t_curr <= 170);
                is_fatigue = (t_curr > 170 && t_curr <= 210);
            elseif k == 2
                is_creep = (t_curr > 312 && t_curr <= 372);
                is_fatigue = (t_curr > 372 && t_curr <= 412);
            elseif k == 3
                is_creep = (t_curr > 20 && t_curr <= 80);
                is_fatigue = (t_curr > 80 && t_curr <= 120);
            elseif k == 4
                is_creep = (t_curr > 25 && t_curr <= 85);
                is_fatigue = (t_curr > 85 && t_curr <= 125);
            end
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + s_dot*dt;
            end
        end
        Df = (Dfc/log(1-Sc/Sg))*log(1-min(Delta_Sf, Sg-1e-5)/Sg);
        pred_lives(k) = 1/(Dc + Df);
    end
end

%% 4. PROFESSIONAL JOURNAL PALETTE
journal_colors = {[0.15, 0.25, 0.55], [0.65, 0.05, 0.05], [0.05, 0.40, 0.20], [0.85, 0.55, 0.10]};
mkr = {'o', 's', '^', 'd'}; styles = {'-', '--', ':', '-.'};

%% 5. PLOT (c): LIFE vs MEAN STRAIN (POINTS + FIT)
figure('Color','w','Position', [100 100 850 850]); hold on;

% Quadratic fit for Mean Strain vs log10(Life)
p_mean = polyfit(mean_strains, log10(pred_lives), 2);
x_fit_mean = linspace(min(mean_strains), max(mean_strains), 100);
y_fit_mean = 10.^(polyval(p_mean, x_fit_mean));

% Plot Fit Line
plot(x_fit_mean, y_fit_mean, 'k--', 'LineWidth', 2.0);

% Plot Points
for i = 1:4
    plot(mean_strains(i), pred_lives(i), mkr{i}, 'MarkerSize', 14, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', journal_colors{i}, 'LineWidth', 1.5);
end

set(gca, 'YScale', 'log', 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Mean Strain, $\epsilon_{mean}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid off; box on; axis square;
xlim([-0.002, 0.017]); ylim([min(pred_lives)*0.5, max(pred_lives)*2]);
text(0.02, 0.98, '(c)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%% 6. PLOT (d): LIFE vs R-RATIO (POINTS + FIT)
figure('Color','w','Position', [100 100 850 850]); hold on;

% Quadratic fit for R-ratio vs log10(Life)
p_r = polyfit(r_vals, log10(pred_lives), 2);
x_fit_r = linspace(min(r_vals), max(r_vals), 100);
y_fit_r = 10.^(polyval(p_r, x_fit_r));

% Plot Fit Line
plot(x_fit_r, y_fit_r, 'k--', 'LineWidth', 2.0);

% Plot Points
for i = 1:4
    plot(r_vals(i), pred_lives(i), mkr{i}, 'MarkerSize', 14, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', journal_colors{i}, 'LineWidth', 1.5);
end

set(gca, 'YScale', 'log', 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid off; box on; axis square;
xlim([-1.2, 0.4]); ylim([min(pred_lives)*0.5, max(pred_lives)*2]);
text(0.02, 0.98, '(d)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold off;
figure('Color','w','Position', [100 100 850 850]); hold on;
bd = bar(1:4, pred_lives, bar_width, 'EdgeColor', 'k', 'LineWidth', 1.5, 'BaseValue', 1);
bd.FaceColor = 'flat';
for i = 1:4, bd.CData(i,:) = journal_colors{i}; end
for i = 1:4
    text(i, pred_lives(i)*1.15, sprintf('%.1f', pred_lives(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 20, 'FontName', 'Times New Roman');
end
set(gca, 'YScale', 'log', 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid on; box on; axis square; 
set(gca, 'XTick', 1:4, 'XTickLabel', r_labels);
xlim([0.5, 4.5]); ylim([1, max(pred_lives)*20]);
text(0.02, 0.98, '(d)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hold off;
%% 5. PLOT (c): LIFE vs R-RATIO (POINTS + FIT)
% ============================================================
figure('Color','w','Position', [100 100 850 850]); hold on;
% Quadratic fit for R-ratio vs log10(Life)
p_r_fit = polyfit(r_vals, log10(pred_lives), 2);
x_range = linspace(min(r_vals), max(r_vals), 100);
y_range = 10.^(polyval(p_r_fit, x_range));

% Plot Fit Line
plot(x_range, y_range, 'k--', 'LineWidth', 2.0, 'HandleVisibility', 'off');

% Plot Points with Legend Entries
h_pts = [];
for i = 1:4
    h_pts(i) = plot(r_vals(i), pred_lives(i), mkr{i}, 'MarkerSize', 14, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', journal_colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('R_{\\epsilon} = %s', r_labels{i}));
end

set(gca, 'YScale', 'log', 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid off; box on; axis square;
xlim([-1.2, 0.4]); ylim([min(pred_lives)*0.5, max(pred_lives)*2]);

% --- ADD LEGEND WITH HOLD TIME ---
lgd = legend(h_pts, 'Location', 'northeast', 'FontSize', 18);
title(lgd, 'Hold Time: 60s'); % Explicitly shows the 60s hold context
legend boxon;

text(0.02, 0.98, '(c)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%% 6. PLOT (d): BAR CHART - LIFE vs R-RATIO
% ============================================================
figure('Color','w','Position', [100 100 850 850]); hold on;
bd = bar(1:4, pred_lives, 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5, 'BaseValue', 1);
bd.FaceColor = 'flat';
for i = 1:4, bd.CData(i,:) = journal_colors{i}; end

% Add values above bars
for i = 1:4
    text(i, pred_lives(i)*1.15, sprintf('%.1f', pred_lives(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 20, 'FontName', 'Times New Roman');
end

set(gca, 'YScale', 'log', 'FontSize', 24, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Strain Ratio, $R_{\epsilon}$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Predicted Life, $N_{f,sim}$', 'Interpreter', 'latex', 'FontSize', 26);
grid on; box on; axis square; 
set(gca, 'XTick', 1:4, 'XTickLabel', r_labels);
xlim([0.5, 4.5]); ylim([1, max(pred_lives)*20]);

% Add text for Hold Time context
text(0.95, 0.05, 't_{hold} = 60 s', 'Units', 'normalized', 'FontSize', 20, ...
    'HorizontalAlignment', 'right', 'FontName', 'Times New Roman');

text(0.02, 0.98, '(d)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold off;