clc; clear; close all;

%% 1. SETUP: DEFINE THE FIVE CASES
cases(1).filename = 'shi_30_30_2_4.csv'; cases(1).epsilon_a = 0.012; cases(1).exp_life = 5;           
cases(2).filename = 'pin_30_30_2_0.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 250;         
cases(3).filename = 'Fatigue.csv';       cases(3).epsilon_a = 0.009; cases(3).exp_life = 370;          
cases(4).filename = 'shi_30_30_1_6.csv'; cases(4).epsilon_a = 0.008; cases(4).exp_life = 562;         
cases(5).filename = '30_30_0_007.csv';   cases(5).epsilon_a = 0.007; cases(5).exp_life = 1455;                 

% GLOBAL CONSTANTS (FIXED)
Q = 6.97e-19; % <--- FIXED as requested
T = 1033; 
R = 8.314; 
D0 = 0; 

%% 2. OPTIMIZATION SETUP
% Parameters x: [n1, Sg, Sc_ratio, Dfc, B1_case1...5]
x0 = [0.6, 0.38, 0.9, 0.9, 3.8e5, 3.8e5, 3.8e5, 3.8e5, 3.8e5]; 

options = optimset('Display', 'iter', 'MaxFunEvals', 8000, 'MaxIter', 8000);
exp_lives = [cases.exp_life];

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

% Cost function now receives Q as a fixed external constant
objFun = @(x) cost_function_n1_variable(x, cases, Q, R, T, D0, exp_lives);

fprintf('\nStarting Optimization (Varying n1 and Damage Parameters)...\n');
[x_opt, fval] = fminsearch(objFun, x0, options);

% Extract Optimized Values
n1_opt = x_opt(1);
Sg_opt = x_opt(2);
Sc_opt = x_opt(3) * Sg_opt;
Dfc_opt = x_opt(4);
B1_opt_vec = x_opt(5:9);

%% 3. FINAL CALCULATION & DAMAGE PARTITIONING
% =========================================================================
pred_lives_opt = zeros(1,5);
Dc_vals = zeros(1,5);
Df_vals = zeros(1,5);

for k = 1:5
    [Dc, Delta_Sf] = calculate_case_damage(cases(k), B1_opt_vec(k), n1_opt, Q, R, T, k);
    
    if Delta_Sf >= Sg_opt
        Df = 1.0;
    else
        den = log(1 - Sc_opt/Sg_opt);
        val_in = 1 - Delta_Sf / Sg_opt;
        log_t = log(max(1e-10, val_in));
        Df = D0 + ((Dfc_opt - D0) / den) * log_t;
    end
    
    Dc_vals(k) = Dc;
    Df_vals(k) = max(0, Df);
    
    total_D = Dc_vals(k) + Df_vals(k);
    pred_lives_opt(k) = 1 / total_D;
end

%% 4. DISPLAY RESULTS
% =========================================================================
fprintf('\n============================================================\n');
fprintf('     FINAL OPTIMIZED PARAMETERS (Q FIXED)                   \n');
fprintf('============================================================\n');
fprintf('Optimized n1:               %.4f\n', n1_opt);
fprintf('Fixed Q:                    %.4e\n', Q);
fprintf('Fatigue Limit (Sg):         %.6f\n', Sg_opt);
fprintf('Critical Damage Stress (Sc): %.6f\n', Sc_opt);
fprintf('Critical Fatigue Damage (Dfc): %.4f\n', Dfc_opt);
fprintf('------------------------------------------------------------\n');

Error_Pct = abs(pred_lives_opt - exp_lives) ./ exp_lives * 100;
T_results = table((1:5)', [cases.epsilon_a]', exp_lives', pred_lives_opt', ...
    Dc_vals', Df_vals', Error_Pct', ...
    'VariableNames', {'Case', 'Strain', 'Exp_Life', 'Pred_Life', 'Creep_Dc', 'Fatigue_Df', 'Error_Pct'});
disp(T_results);

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
        
        % Dwell/Cycle Logic for 30/30 Cases (Same as original)
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

function err = cost_function_n1_variable(x, cases, Q, R, T, D0, exp_lives)
    n1_c = x(1); Sg_c = x(2); Sc_r = x(3); Dfc_c = x(4); B1_v = x(5:9); 
    
    penalty = 0;
    % Constraints
    if n1_c < 0.2 || n1_c > 0.8, penalty = penalty + 1e7; end
    if Sg_c <= 0.001 || Sg_c > 1.5, penalty = penalty + 1e7; end
    if Sc_r <= 0.01 || Sc_r >= 0.99, penalty = penalty + 1e7; end
    if any(B1_v < 1e3), penalty = penalty + 1e7; end
    
    if penalty > 0, err = 1e10 + penalty; return; end
    
    Sc_c = Sc_r * Sg_c;
    p_lives = zeros(1,5);
    for k = 1:5
        [Dc, Delta_Sf] = calculate_case_damage(cases(k), B1_v(k), n1_c, Q, R, T, k);
        if Delta_Sf >= Sg_c, Df = 1.0; else
            den = log(1 - Sc_c/Sg_c);
            val_in = 1 - Delta_Sf / Sg_c;
            Df = D0 + ((Dfc_c - D0) / den) * log(max(1e-10, val_in));
        end
        total_D = Dc + max(0, Df);
        p_lives(k) = 1 / max(1e-12, total_D);
    end
    err = sqrt(mean((log10(exp_lives) - log10(p_lives)).^2)) + penalty;
end