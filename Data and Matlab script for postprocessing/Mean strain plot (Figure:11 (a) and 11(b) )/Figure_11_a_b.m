clc; clear; close all;


%%%%%%% Figure -11%%%%%%%%%% 


%% 1. SETTINGS & FILENAMES
files = {'mean_strain_0.csv', 'mean_strain_0_002.csv', ...
         'mean_strain_0_01.csv', 'mean_strain_0_015.csv'};                    
starts = [1900, 1500, 1001, 973]; 
mean_strains = [0, 0.002, 0.010, 0.015]; 
eps_a_vals = [0.010, 0.002, 0.010, 0.015]; 

% Standardized labels for legends
labels = {'R_{\epsilon} = -1', 'R_{\epsilon} = -0.67', 'R_{\epsilon} = 0', 'R_{\epsilon} = 0.2'};

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


%% 4. FORMATTING (Figure 1: Stress-Strain)
set(gca, 'FontSize', 22, ... % Synced with Figure 2
         'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, ...
         'TickDir', 'in');

xlabel('Strain, $\epsilon$', 'Interpreter','latex', 'FontSize', 26);
ylabel('Stress, $\sigma$ (MPa)', 'Interpreter','latex', 'FontSize', 26);

% Legend moved to SOUTHWEST to match Figure 2
lgd1 = legend(labels, 'Location','southwest');
lgd1.FontSize = 24; 
lgd1.FontName = 'Times New Roman'; 
legend boxon;

% Annotation Box in SOUTHEAST
text(0.0285, -1150, '(60/0) hold ', ...
    'FontSize', 24, ... 
    'FontName', 'Times New Roman', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'HorizontalAlignment', 'right', ... 
    'VerticalAlignment', 'bottom');

xlim([-0.03 0.03]);
ylim([-1200 1200]);
grid on; % Synced with Figure 2
box on;
axis square;
hold off;



%% 7. PLOT (c): STRAIN vs TIME
figure('Color','w','Position', [100 100 850 850]); hold on;
h_time = [];
for i = 1:4
    if ~isempty(processed_data{i})
        % Define data for marker placement
        t_p = processed_data{i}.t_full;
        e_p = processed_data{i}.e_full;
        
        % Plot with markers matching Figure 1 styling
        h_time(i) = plot(t_p, e_p, ...
            'Color', journal_colors{i}, ...
            'LineStyle', styles{i}, ...
            'LineWidth', 2.0, ...
            'Marker', markers{i}, ...
            'MarkerSize', 9, ...
            'MarkerFaceColor', 'w', ...
            'MarkerIndices', round(linspace(1, length(t_p), 30))); % 30 markers for clarity
    end
end

set(gca, 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickDir', 'in');
xlabel('Time, $t$ (s)', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Strain, $\epsilon$', 'Interpreter', 'latex', 'FontSize', 26);
grid off; box on; axis square; 

% Set limits
xlim([0 210]); 
ylim([-0.025 0.03]); 

% --- Annotation Box in SOUTHEAST ---
text(205, -0.0236, '(60/0) hold ', ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'HorizontalAlignment', 'right', ... 
    'VerticalAlignment', 'bottom');

% Legend in SOUTHWEST
lgd3 = legend(h_time, labels, 'Location', 'southwest');
lgd3.FontSize = 24; lgd3.FontName = 'Times New Roman'; legend boxon;
hold off;

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

%% 9. SAVE FIGURES AS PDFs
% Save Figure 1 (Stress-Strain Hysteresis)
fig1 = figure(1); 
exportgraphics(fig1, 'Stress_Strain_Hysteresis.pdf', 'ContentType', 'vector');

% Save Figure 2 (Strain vs Time)
fig2 = figure(2);
exportgraphics(fig2, 'Strain_vs_Time.pdf', 'ContentType', 'vector');

fprintf('Figures saved successfully as PDFs.\n');
