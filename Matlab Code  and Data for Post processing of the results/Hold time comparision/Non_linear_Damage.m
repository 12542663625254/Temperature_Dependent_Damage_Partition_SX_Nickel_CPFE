clc; clear; close all;

%% 1. SETUP: DEFINE MATERIAL CONSTANTS (980°C)
Sg = 0.055; Sc = 0.01268; Dfc = 0.055;
n1 = 0.6; Q = 6.97e-19; T = 1253; R = 8.314;
q_nlds = 0.576; % Exponent for Non-linear Damage Summation

%% 2. DEFINE DATASETS
ds(1).label = '120/0 Dwell';
ds(1).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(1).files = {'120_0_2_4.csv', '120_0_2_0.csv', '120_0_1_8.csv', '120_0_1_6.csv'};

ds(2).label = '60/0 Dwell';
ds(2).epsilon_a = [0.010, 0.010, 0.010, 0.010];
ds(2).files = {'Fatigue_60_0_2_0.csv', 'mean_strain_0_002.csv', 'mean_strain_0_01.csv', 'mean_strain_0_015.csv'};

ds(3).label = '240/0 Dwell';
ds(3).epsilon_a = [0.012, 0.010, 0.009, 0.008];
ds(3).files = {'240_0_2_4.csv', '240_0_2_0.csv', '240_0_1_8.csv', '240_0_1_6.csv'};

% Initialize cell arrays to store results
all_pred_lives_LDS = cell(1, length(ds));
all_pred_lives_NLDS = cell(1, length(ds));

%% 3. MAIN PROCESSING LOOP
for d = 1:length(ds)
    num_cases = length(ds(d).epsilon_a);
    pred_lives_lds = zeros(1, num_cases);
    pred_lives_nlds = zeros(1, num_cases);
    
    for k = 1:num_cases
        eps_a = ds(d).epsilon_a(k);
        current_file = ds(d).files{k};
        B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
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
                
                % Time Logic Separation
                if d == 1 % 120/0
                    if k==1, if (t_curr > 12 && t_curr <= 132), is_creep = true; elseif (t_curr > 132 && t_curr <= 180), is_fatigue = true; end
                    elseif k==2, if (t_curr > 170 && t_curr <= 290), is_creep = true; elseif (t_curr > 290 && t_curr <= 330), is_fatigue = true; end
                    elseif k==3, if (t_curr > 9 && t_curr <= 129), is_creep = true; elseif (t_curr > 129 && t_curr <= 165), is_fatigue = true; end
                    elseif k==4, if (t_curr > 8 && t_curr <= 128), is_creep = true; elseif (t_curr > 128 && t_curr <= 160), is_fatigue = true; end
                    end
                elseif d == 2 % 60/0
                    if k==1, if (t_curr > 10 && t_curr <= 70), is_creep = true; elseif (t_curr > 70 && t_curr <= 110), is_fatigue = true; end
                    elseif k==2, if (t_curr > 12 && t_curr <= 72), is_creep = true; elseif (t_curr > 72 && t_curr <= 112), is_fatigue = true; end
                    elseif k==3, if (t_curr > 20 && t_curr <= 80), is_creep = true; elseif (t_curr > 80 && t_curr <= 120), is_fatigue = true; end
                    elseif k==4, if (t_curr > 25 && t_curr <= 85), is_creep = true; elseif (t_curr > 85 && t_curr <= 125), is_fatigue = true; end
                    end
                elseif d == 3 % 240/0
                    if k==1, if (t_curr > 300 && t_curr <= 540), is_creep = true; elseif (t_curr > 540 && t_curr <= 828), is_fatigue = true; end
                    elseif k==2, if (t_curr > 10 && t_curr <= 250), is_creep = true; elseif (t_curr > 250 && t_curr <= 290), is_fatigue = true; end
                    elseif k==3, if (t_curr > 285 && t_curr <= 525), is_creep = true; elseif (t_curr > 525 && t_curr <= 561), is_fatigue = true; end
                    elseif k==4, if (t_curr > 280 && t_curr <= 520), is_creep = true; elseif (t_curr > 520 && t_curr <= 552), is_fatigue = true; end
                    end
                end
                
                if is_creep && s_dot > 0, Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
                elseif is_fatigue, Delta_Sf = Delta_Sf + s_dot*dt; end
            end
            Df = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
            
            % --- CALCULATE LIVES ---
            pred_lives_lds(k) = 1 / (Dc + Df);                 % Linear Rule
            pred_lives_nlds(k) = 1 / (Dc^q_nlds + Df^q_nlds); % Non-linear Rule
        else
            pred_lives_lds(k) = NaN; pred_lives_nlds(k) = NaN;
        end
    end
    all_pred_lives_LDS{d} = pred_lives_lds;
    all_pred_lives_NLDS{d} = pred_lives_nlds;
end

%% 4. FIGURE: CYCLES TO FAILURE vs HOLD TIME (NON-LINEAR PLOT)
hold_times_labels = {'60s', '120s', '240s'};
strain_legend = {'$\Delta\epsilon = 2.4\%$', '$\Delta\epsilon = 2.0\%$', '$\Delta\epsilon = 1.8\%$', '$\Delta\epsilon = 1.6\%$'};

% Life matrix: Rows = Hold Times (60, 120, 240)
life_matrix_nlds = [all_pred_lives_NLDS{2}; all_pred_lives_NLDS{1}; all_pred_lives_NLDS{3}];

figure('Name', 'Life vs Hold Time - NLDS', 'Color','w', 'Position', [100 100 800 800]); 
b_hold = bar(life_matrix_nlds, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.2);

set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'TickDir', 'in'); 
axis square; grid on; box on;

set(gca, 'XTickLabel', hold_times_labels);
xlabel('Hold Duration (seconds)', 'FontName', 'Times New Roman', 'FontSize', 26);
ylabel('Predicted Life, $N_{f}$ (NLDS)', 'Interpreter', 'latex', 'FontSize', 26);
ylim([10 10000]); 

legend(strain_legend, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);

% Numerical labels on bars
for i = 1:size(life_matrix_nlds, 2)
    xtips = b_hold(i).XEndPoints;
    ytips = b_hold(i).YData;
    labels = string(round(ytips));
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontName', 'Times New Roman');
end

% --- TRANSPARENT ANNOTATION BOX FOR NLDS ---
bX = 0.58; bY = 0.65; bW = 0.30; bH = 0.12;
annotation('rectangle', [bX, bY, bW, bH], 'FaceColor', 'none', 'LineWidth', 1.5);
annotation('textbox', [bX, bY, bW, bH], ...
    'String', sprintf('NLDS Rule\n$N_f = 1 / (D_c^q + D_f^q)$\n$q = 0.576$'), ...
    'Interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 18, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');