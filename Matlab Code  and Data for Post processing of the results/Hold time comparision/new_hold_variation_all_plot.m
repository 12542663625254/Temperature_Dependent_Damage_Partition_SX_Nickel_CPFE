clc; clear; close all;

%% 1. SETUP: DEFINE MATERIAL CONSTANTS (980°C)
Sg = 0.38; Sc = 0.9*Sg; Dfc = 0.9;
n1 = 0.6; Q = 6.97e-19; T = 1033; R = 8.314;

%% 2. DEFINE DATASETS
ds(1).label = '120/0 Dwell'; ds(1).files = {'120_0_2_4.csv', '120_0_2_0.csv', '120_0_1_8.csv', '120_0_1_6.csv'};
ds(2).label = '60/0 Dwell';  ds(2).files = {'Fatigue_60_0_2_4.csv', 'Fatigue_60_0_2_0.csv', 'Fatigue_60_0_1_8.csv', 'Fatigue_60_0_0_8_perc.csv'};
ds(3).label = '240/0 Dwell'; ds(3).files = {'240_0_2_4.csv', '240_0_2_0.csv', '240_0_1_8.csv', '240_0_1_6.csv'};
ds(4).label = '300/0 Dwell'; ds(4).files = {'300_0_2_4_new.csv', '300_0_2_0_new.csv', '300_0_1_8_new.csv', '300_0_1_6_new.csv'};

% --- DATASET 5: 0s Hold (Pure Fatigue) ---
ds(5).label = '0s Hold';
ds(5).files = {'cyclic_2_4.csv', 'cyclic_2_0.csv', 'cyclic_1_8.csv', 'cyclic_1_6.csv'}; 

for d = 1:5, ds(d).epsilon_a = [0.012, 0.010, 0.009, 0.008]; end

all_pred_lives = cell(1, 5); 
all_Dc = cell(1, 5); 
all_Df = cell(1, 5);

%% 3. MAIN PROCESSING LOOP
for d = 1:length(ds)
    num_cases = length(ds(d).epsilon_a);
    pred_lives = zeros(1, num_cases);
    dc_vals = zeros(1, num_cases); 
    df_vals = zeros(1, num_cases);
    
    for k = 1:num_cases
        eps_a = ds(d).epsilon_a(k);
        current_file = ds(d).files{k};
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
                
                % Time Logic
                if d == 1 % 120/0
                    if k == 1, if (t_curr > 12 && t_curr <= 132), is_creep = true; elseif (t_curr > 132 && t_curr <= 180), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 170 && t_curr <= 290), is_creep = true; elseif (t_curr > 290 && t_curr <= 330), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 129), is_creep = true; elseif (t_curr > 129 && t_curr <= 165), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 128), is_creep = true; elseif (t_curr > 128 && t_curr <= 160), is_fatigue = true; end
                    end
                elseif d == 2 % 60/0
                    if k == 1, if (t_curr > 12 && t_curr <= 72), is_creep = true; elseif (t_curr > 72 && t_curr <= 120), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 70), is_creep = true; elseif (t_curr > 70 && t_curr <= 110), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 69), is_creep = true; elseif (t_curr > 69 && t_curr <= 105), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 68), is_creep = true; elseif (t_curr > 68 && t_curr <= 100), is_fatigue = true; end
                    end
                elseif d == 3 % 240/0
                    if k == 1, if (t_curr > 300 && t_curr <= 540), is_creep = true; elseif (t_curr > 540 && t_curr <= 828), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 250), is_creep = true; elseif (t_curr > 250 && t_curr <= 290), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 285 && t_curr <= 525), is_creep = true; elseif (t_curr > 525 && t_curr <= 561), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 280 && t_curr <= 520), is_creep = true; elseif (t_curr > 520 && t_curr <= 552), is_fatigue = true; end
                    end
                elseif d == 4 % 300/0
                    if k == 1, if (t_curr > 12 && t_curr <= 312), is_creep = true; elseif (t_curr > 312 && t_curr <= 360), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 310), is_creep = true; elseif (t_curr > 310 && t_curr <= 350), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 309), is_creep = true; elseif (t_curr > 309 && t_curr <= 345), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 308), is_creep = true; elseif (t_curr > 308 && t_curr <= 340), is_fatigue = true; end
                    end
                elseif d == 5 % 0s Hold (Pure Fatigue)
                    if k == 1, if (t_curr > 12 && t_curr <= 60), is_fatigue = true; end
                    elseif k == 2, if (t_curr > 10 && t_curr <= 50), is_fatigue = true; end
                    elseif k == 3, if (t_curr > 9 && t_curr <= 45), is_fatigue = true; end
                    elseif k == 4, if (t_curr > 8 && t_curr <= 40), is_fatigue = true; end
                    end
                end
                
                if is_creep && s_dot > 0, Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
                elseif is_fatigue, Delta_Sf = Delta_Sf + s_dot*dt; end
            end
            
            % --- DAMAGE ASSIGNMENT ---
            if d == 5
                df_vals(k) = 1/500; % Anchored at 500 cycles
                dc_vals(k) = NaN;   % Red line won't plot at 0s
                pred_lives(k) = 500;
            else
                df_vals(k) = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
                dc_vals(k) = Dc;
                pred_lives(k) = 1/(dc_vals(k) + df_vals(k));
            end
        else
            pred_lives(k) = NaN; dc_vals(k) = NaN; df_vals(k) = NaN;
        end
    end
    all_pred_lives{d} = pred_lives;
    all_Dc{d} = dc_vals;
    all_Df{d} = df_vals;
end

%% 4. VISUALIZATION: COMPARISON OF 0.8%, 1.0%, and 1.2% (Hold Periods Only)
hold_durations = [60, 120, 240, 300]; 
idx_mapping = [2, 1, 3, 4]; % Maps 60s, 120s, 240s, 300s holds
num_holds = length(hold_durations);

Dc_08 = zeros(1, num_holds); Df_08 = zeros(1, num_holds);
Dc_10 = zeros(1, num_holds); Df_10 = zeros(1, num_holds);
Dc_12 = zeros(1, num_holds); Df_12 = zeros(1, num_holds);

for i = 1:num_holds
    idx = idx_mapping(i);
    Dc_08(i) = all_Dc{idx}(4); Df_08(i) = all_Df{idx}(4);
    Dc_10(i) = all_Dc{idx}(2); Df_10(i) = all_Df{idx}(2);
    Dc_12(i) = all_Dc{idx}(1); Df_12(i) = all_Df{idx}(1);
end

% Define requested colors
color12 = [1, 0, 0];       % Red
color10 = [0, 0, 1];       % Blue
color08 = [0, 0.5, 0];     % Dark Green

figure('Color','w','Position', [100 100 1200 900]);
hold on; grid on;
y_ticks = [0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1];

% --- Left Axis (Creep Damage - Dashed) ---
yyaxis left
c1 = plot(hold_durations, Dc_12, '-s', 'Color', color12, 'LineWidth', 3, 'MarkerSize', 12, 'MarkerFaceColor', color12);
c2 = plot(hold_durations, Dc_10, '-p', 'Color', color10, 'LineWidth', 3, 'MarkerSize', 12, 'MarkerFaceColor', color10);
c3 = plot(hold_durations, Dc_08, '-o', 'Color', color08, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color08);
ylabel('Creep damage per cycle, $D_c$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'YScale', 'log', 'YColor', 'k', 'YLim', [0.0002, 4e-1]); % Y-axis set to black for neutrality
yticks(y_ticks);

% --- Right Axis (Fatigue Damage - Solid) ---
yyaxis right
f1 = plot(hold_durations, Df_12, '--s', 'Color', color12, 'LineWidth', 3, 'MarkerSize', 12);
f2 = plot(hold_durations, Df_10, '--p', 'Color', color10, 'LineWidth', 3, 'MarkerSize', 12);
f3 = plot(hold_durations, Df_08, '--o', 'Color', color08, 'LineWidth', 3, 'MarkerSize', 10);
ylabel('Fatigue damage per cycle, $D_f$', 'Interpreter', 'latex', 'FontSize', 28);
set(gca, 'YScale', 'log', 'YColor', 'k', 'YLim', [0.0002, 4e-1]); % Y-axis set to black for neutrality
yticks(y_ticks);

% Labels and Formatting
xlim([50 310]);
xlabel('Tensile Hold Time (s)', 'FontSize', 28, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'in', 'Box', 'on');

% Updated Legend with explicit style labeling
lgd = legend([c1, f1, c2, f2, c3, f3], ...
    {'$D_c$ (1.2\%)', ' $D_f$ (1.2\%)', ...
     '$D_c$ (1.0\%)', '$D_f$ (1.0\%)', ...
     '$D_c$ (0.8\%)', '$D_f$ (0.8\%)'}, ...
    'Location', 'northoutside', 'Interpreter', 'latex', 'FontSize', 18, 'NumColumns', 3);

axis square;
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off');