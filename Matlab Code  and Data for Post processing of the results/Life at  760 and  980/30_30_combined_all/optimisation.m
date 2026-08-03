clc; clear; close all;

%% 1. SETUP: DEFINE THE FIVE CASES (30/30 Dwell)
% =========================================================================
% Case data for DD6 Nickel-Based Superalloy at 760°C
cases(1).filename = 'shi_30_30_2_4.csv'; cases(1).epsilon_a = 0.012; cases(1).exp_life = 5;           
cases(2).filename = 'pin_30_30_2_0.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 250;         
cases(3).filename = 'Fatigue.csv';       cases(3).epsilon_a = 0.009; cases(3).exp_life = 370;          
cases(4).filename = 'shi_30_30_1_6.csv'; cases(4).epsilon_a = 0.008; cases(4).exp_life = 562;         
cases(5).filename = '30_30_0_007.csv';   cases(5).epsilon_a = 0.007; cases(5).exp_life = 1455;                 

% GLOBAL CONSTANTS
n1 = 0.6; 
Q = 6.97e-19; 
T = 1033; 
R = 8.314; 
D0 = 0; 

%% 2. OPTIMIZATION SETUP
% =========================================================================
% Parameters x: [Sg, Sc_ratio, Dfc, B1_case1, B1_case2, B1_case3, B1_case4, B1_case5]
x0 = [0.38, 0.9, 0.9, 3.8e5, 3.8e5, 3.8e5, 3.8e5, 3.8e5]; 
options = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'MaxIter', 5000);
exp_lives = [cases.exp_life];

% Pre-load data to speed up optimization
fprintf('Pre-loading CSV data...\n');
for k = 1:5
    if isfile(cases(k).filename)
        raw = readmatrix(cases(k).filename);
        cases(k).time = raw(:, 1);
        cases(k).s_dot = raw(:, 25);
    else
        error('File %s not found!', cases(k).filename);
    end
end

objFun = @(x) cost_function_B_variable(x, cases, n1, Q, R, T, D0, exp_lives);
fprintf('\nStarting Optimization with Strain-Specific B1 values...\n');
[x_opt, fval] = fminsearch(objFun, x0, options);

% Extract Optimized Values
Sg_opt = x_opt(1);
Sc_ratio_opt = x_opt(2);
Sc_opt = Sc_ratio_opt * Sg_opt;
Dfc_opt = x_opt(3);
B1_opt_vec = x_opt(4:8);

%% 3. FINAL CALCULATION (Damage Partitioning)
% =========================================================================
pred_lives_opt = zeros(1,5);
Dc_vals = zeros(1,5);
Df_vals = zeros(1,5);

for k = 1:5
    [Dc, Delta_Sf] = calculate_case_damage(cases(k), B1_opt_vec(k), n1, Q, R, T, k);
    
    % Fatigue Damage Calculation based on Critical Stress Sg
    if Delta_Sf >= Sg_opt
        Df = 1.0;
    else
        den = log(1 - Sc_opt/Sg_opt);
        val_in = 1 - Delta_Sf / Sg_opt;
        
        if val_in <= 1e-9
            log_t = -50;
        else
            log_t = log(max(1e-10, val_in));
        end
        
        Df = D0 + ((Dfc_opt - D0) / den) * log_t;
    end
    
    Df = max(0, Df);
    Dc_vals(k) = Dc;
    Df_vals(k) = Df;
    
    total_D = Dc + Df;
    
    if total_D <= 1e-9
        pred_lives_opt(k) = 1e9;
    else
        pred_lives_opt(k) = 1 / total_D;
    end
end

%% 4. DISPLAY RESULTS & DAMAGE ANALYSIS
% =========================================================================
fprintf('\n============================================================\n');
fprintf('     FINAL OPTIMIZED PARAMETERS & DAMAGE RESULTS            \n');
fprintf('============================================================\n');
fprintf('Fatigue Limit Parameter (Sg):       %.6f\n', Sg_opt);
fprintf('Critical Damage Stress Ratio (Sc):  %.6f\n', Sc_opt);
fprintf('Critical Fatigue Damage (Dfc):      %.4f\n', Dfc_opt);
fprintf('------------------------------------------------------------\n');

% Create detailed comparison table
Error_Pct = abs(pred_lives_opt - exp_lives) ./ exp_lives * 100;
T_results = table((1:5)', [cases.epsilon_a]', exp_lives', pred_lives_opt', ...
    Dc_vals', Df_vals', Error_Pct', ...
    'VariableNames', {'Case', 'Strain_Amp', 'Exp_Life', 'Pred_Life', 'Creep_Dc', 'Fatigue_Df', 'Error_Pct'});

disp(T_results);

fprintf('\nRelative Damage Contribution per Cycle:\n');
for k = 1:5
    total = Dc_vals(k) + Df_vals(k);
    fprintf('  Case %d (ea=%.3f): Creep = %.2f%%, Fatigue = %.2f%%\n', ...
        k, cases(k).epsilon_a, (Dc_vals(k)/total)*100, (Df_vals(k)/total)*100);
end

%% 5. PLOT RESULTS
% =========================================================================
figure('Color', 'w', 'Position', [100, 100, 900, 700]);
all_vals = [exp_lives, pred_lives_opt];
limits = [min(all_vals(all_vals>0))*0.5, max(all_vals)*2];
loglog(limits, limits, 'k-', 'LineWidth', 2); hold on;
loglog(limits, limits*2, 'k--', 'LineWidth', 1.5);
loglog(limits, limits*0.5, 'k--', 'LineWidth', 1.5);

colors = {'r', 'b', 'g', 'm', 'k'}; 
for k = 1:5
    loglog(exp_lives(k), pred_lives_opt(k), 's', ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k}, ...
        'DisplayName', sprintf('Strain %.3f', cases(k).epsilon_a));
end
grid on; axis square; xlim(limits); ylim(limits);
xlabel('Experimental Life (Cycles)'); ylabel('Predicted Life (Cycles)');
title('Life Prediction with Creep-Fatigue Partitioning');
legend('Location', 'northwest');

%% HELPER FUNCTIONS
function [Dc_accum, Delta_Sf_accum] = calculate_case_damage(case_struct, B1, n1, Q, R, T, k)
    phi = B1 * exp(-Q / (R * T));
    t = case_struct.time;
    s_dot_vec = case_struct.s_dot;
    Dc_accum = 0; Delta_Sf_accum = 0;
    
    for i = 2:length(t)
        dt = t(i) - t(i-1);
        s_dot = s_dot_vec(i);
        t_curr = t(i);
        is_creep = false; is_fatigue = false;
        
        % Dwell/Cycle Logic for 30/30 Cases
        if k == 1
            if (t_curr>12 && t_curr<=42)||(t_curr>66 && t_curr<=96), is_creep = true;
            elseif (t_curr>42 && t_curr<=66)||(t_curr>96 && t_curr<=120), is_fatigue = true; end
        elseif k == 2
            if (t_curr>10 && t_curr<=40)||(t_curr>60 && t_curr<=90), is_creep = true;
            elseif (t_curr>40 && t_curr<=60)||(t_curr>90 && t_curr<=110), is_fatigue = true; end
        elseif k == 3
            if (t_curr>681 && t_curr<=711)||(t_curr>729 && t_curr<=759), is_creep = true;
            elseif (t_curr>711 && t_curr<=729)||(t_curr>759 && t_curr<=777), is_fatigue = true; end
        elseif k == 4
            if (t_curr>560 && t_curr<=590)||(t_curr>606 && t_curr<=636), is_creep = true;
            elseif (t_curr>590 && t_curr<=606)||(t_curr>636 && t_curr<=652), is_fatigue = true; end
        elseif k == 5
            if (t_curr>7 && t_curr<=37)||(t_curr>51 && t_curr<=81), is_creep = true;
            elseif (t_curr>37 && t_curr<=51)||(t_curr>81 && t_curr<=95), is_fatigue = true; end
        end
        
        if is_creep && s_dot > 0
            Dc_accum = Dc_accum + ((1/phi) * (s_dot)^(1-n1) * dt);
        elseif is_fatigue
            Delta_Sf_accum = Delta_Sf_accum + (s_dot * dt);
        end
    end
end

function err = cost_function_B_variable(x, cases, n1, Q, R, T, D0, exp_lives)
    Sg_curr = x(1); Sc_ratio = x(2); Dfc_curr = x(3); B1_vec = x(4:8); 
    penalty = 0;
    if Sg_curr <= 0.001 || Sg_curr > 1.0, penalty = penalty + 1e7; end
    if Sc_ratio <= 0.01 || Sc_ratio >= 0.999, penalty = penalty + 1e7; end
    if Dfc_curr <= 0.001 || Dfc_curr > 10, penalty = penalty + 1e7; end
    if any(B1_vec < 1e2) || any(B1_vec > 1e9), penalty = penalty + 1e7; end
    if penalty > 0, err = 1e10 + penalty; return; end
    
    Sc_curr = Sc_ratio * Sg_curr;
    p_lives = zeros(1,5);
    for k = 1:5
        [Dc, Delta_Sf] = calculate_case_damage(cases(k), B1_vec(k), n1, Q, R, T, k);
        if Delta_Sf >= Sg_curr
            Df = 1.0;
        else
            den = log(1 - Sc_curr/Sg_curr);
            val_in = 1 - Delta_Sf / Sg_curr;
            if val_in <= 1e-9, log_t = -50; else, log_t = log(max(1e-10, val_in)); end
            Df = D0 + ((Dfc_curr - D0) / den) * log_t;
        end
        total_D = Dc + max(0, Df);
        if total_D <= 1e-9, p_lives(k) = 1e9; else, p_lives(k) = 1 / total_D; end
    end
    err = sqrt(mean((log10(exp_lives) - log10(p_lives)).^2)) + penalty;
end