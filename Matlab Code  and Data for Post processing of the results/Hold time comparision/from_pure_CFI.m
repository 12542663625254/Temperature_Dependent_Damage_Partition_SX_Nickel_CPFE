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



%% 6. NEW FIGURE: MODIFIED LIFE vs HOLD TIME (Including 0s Hold)
% ========================================================================
hold_times_labels_mod = {'0s', '60s', '120s', '240s', '300s'};
strain_legend_mod = {'$\Delta\epsilon = 2.4\%$', '$\Delta\epsilon = 2.0\%$', '$\Delta\epsilon = 1.8\%$', '$\Delta\epsilon = 1.6\%$'};

% --- DATA PREPARATION ---
life_0s = [55, 222, 612, 2527]; 
life_matrix_existing = [
    all_pred_lives{2}; ... % 60s
    all_pred_lives{1}; ... % 120s
    all_pred_lives{3}; ... % 240s
    all_pred_lives{4}      % 300s
];

life_matrix_existing(1, 2) = 100; 
life_matrix_final = [life_0s; life_matrix_existing];

% --- PLOTTING ---
% Increased width (1100) to make the plot "wider" overall
figure('Name', 'Life vs Hold Time with 0s', 'Color','w', 'Position', [100 150 1100 800]); 

% Changed width to 0.85 for wider individual bars
b_mod = bar(life_matrix_final, 0.85, 'grouped', 'EdgeColor', 'k', 'LineWidth', 1.2);

% --- STYLE FORMATTING ---
set(gca, 'YScale', 'log', 'FontSize', 22, 'FontName', 'Times New Roman', 'TickDir', 'in');
% Removed 'axis square' to allow the wide aspect ratio to take effect
grid off; 
box on;
set(gca, 'XTickLabel', hold_times_labels_mod);
xlabel('Tensile Holding Time (s)', 'FontName', 'Times New Roman', 'FontSize', 26);
ylabel('Cycles to Failure ($N_{f}$)', 'Interpreter', 'latex', 'FontSize', 26);

ylim([10 10000]); 
legend(strain_legend_mod, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 16);

% --- NUMERICAL LABELS ON BARS ---
for i = 1:size(life_matrix_final, 2)
    xtips = b_mod(i).XEndPoints;
    ytips = b_mod(i).YData;
    labels = string(round(ytips));
    
    % Adjusted font size and alignment for wider bars
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 25, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end
axis square;
text(0.02, 0.95, '(d)', 'Units', 'normalized', 'FontSize', 30, ...
    'FontName', 'Times New Roman', 'FontWeight', 'bold');
hold off;