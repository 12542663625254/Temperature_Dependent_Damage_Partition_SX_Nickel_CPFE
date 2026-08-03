clc; clear; close all;

%% 1. SETUP: DEFINE THE FIVE CASES
% =========================================================================
cases(1).filename = 'shi_30_30_2_4.csv'; cases(1).epsilon_a = 0.012; cases(1).exp_life = 5;           
cases(2).filename = 'pin_30_30_2_0.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 250;         
cases(3).filename = 'Fatigue.csv';       cases(3).epsilon_a = 0.009; cases(3).exp_life = 370;          
cases(4).filename = 'shi_30_30_1_6.csv'; cases(4).epsilon_a = 0.008; cases(4).exp_life = 562;         
cases(5).filename = '30_30_0_007.csv';   cases(5).epsilon_a = 0.007; cases(5).exp_life = 1455;                 

% GLOBAL CONSTANTS
n1 = 0.6; Q = 6.97e-19; T = 1033; R = 8.314; D0 = 0; 
Sc_ratio_fixed = 0.9; Dfc_fixed = 0.9;

%% 2. OPTIMIZATION SETUP
% =========================================================================
x0 = [0.8, 1e7, 1e7, 1e7, 1e7, 1e7]; 
options = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'MaxIter', 5000);

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

fprintf('\nOptimizing: Matching Life while forcing Df ≈ 10*Dc...\n');
[x_opt, fval] = fminsearch(@(x) cost_function_entropy_scaled(x, cases, n1, Q, R, T, Sc_ratio_fixed, Dfc_fixed), x0, options);

Sg_opt = x_opt(1);
B1_opt_vec = x_opt(2:end);

%% 3. FINAL EVALUATION & ENTROPY COMPARISON
% =========================================================================
pred_lives = zeros(1,5);
Dc_vals = zeros(1,5); Df_vals = zeros(1,5);
S_creep = zeros(1,5); S_fatigue = zeros(1,5);

for k = 1:5
    [Dc, Delta_Sf, Sc_val, Sf_val] = calculate_advanced_metrics(cases(k), B1_opt_vec(k), n1, Q, R, T, k);
    
    den = log(1 - Sc_ratio_fixed); 
    val_in = 1 - Delta_Sf / Sg_opt;
    log_t = log(max(1e-10, val_in));
    Df = max(1e-12, D0 + ((Dfc_fixed - D0) / den) * log_t);
    
    Dc_vals(k) = Dc; Df_vals(k) = Df;
    S_creep(k) = Sc_val; S_fatigue(k) = Sf_val;
    pred_lives(k) = 1 / (Dc + Df);
end

%% 4. DISPLAY RESULTS
% =========================================================================
fprintf('\n============================================================\n');
fprintf('     THERMODYNAMIC & DAMAGE PARTITIONING RESULTS            \n');
fprintf('============================================================\n');

Ratio_Df_Dc = Df_vals ./ Dc_vals;
T_results = table((1:5)', [cases.epsilon_a]', Dc_vals', Df_vals', Ratio_Df_Dc', ...
    S_creep', S_fatigue', pred_lives', ...
    'VariableNames', {'Case','Strain','Dc','Df','Ratio_Df_Dc','Sgen_Creep','Sgen_Fat','Pred_Life'});
disp(T_results);

%% HELPER FUNCTIONS
function [Dc, Delta_Sf, S_creep, S_fatigue] = calculate_advanced_metrics(case_struct, B1, n1, Q, R, T, k)
    phi = B1 * exp(-Q / (R * T));
    t = case_struct.time; s_dot = case_struct.s_dot;
    Dc = 0; Delta_Sf = 0; S_creep = 0; S_fatigue = 0;
    
    for i = 2:length(t)
        dt = t(i) - t(i-1);
        ds_dot = s_dot(i);
        t_curr = t(i);
        is_c = false; is_f = false;
        
        % Full Dwell/Cycle Logic for all 5 Cases
        if k == 1
            if (t_curr>228 && t_curr<=258)||(t_curr>282&& t_curr<=312), is_c = true;
            elseif (t_curr>258 && t_curr<=282)||(t_curr>312 && t_curr<=336), is_f = true; end
        elseif k == 2
            if (t_curr>10 && t_curr<=40)||(t_curr>60 && t_curr<=90), is_c = true;
            elseif (t_curr>40 && t_curr<=60)||(t_curr>90 && t_curr<=110), is_f = true; end
        elseif k == 3
            if (t_curr>9 && t_curr<=39)||(t_curr>57 && t_curr<=87), is_c = true;
            elseif (t_curr>39 && t_curr<=57)||(t_curr>87 && t_curr<=105), is_f = true; end
        elseif k == 4
            if (t_curr>560 && t_curr<=590)||(t_curr>606 && t_curr<=636), is_c = true;
            elseif (t_curr>590 && t_curr<=606)||(t_curr>636 && t_curr<=652), is_f = true; end
        elseif k == 5
            if (t_curr>7 && t_curr<=37)||(t_curr>51 && t_curr<=81), is_c = true;
            elseif (t_curr>37 && t_curr<=51)||(t_curr>81 && t_curr<=95), is_f = true; end
        end
        
        dissipation = abs(ds_dot); % Simplified Thermodynamic dissipation power
        if is_c && ds_dot > 0
            Dc = Dc + ((1/phi) * (ds_dot)^(1-n1) * dt);
            S_creep = S_creep + (dissipation / T) * dt;
        elseif is_f
            Delta_Sf = Delta_Sf + (ds_dot * dt);
            S_fatigue = S_fatigue + (dissipation / T) * dt;
        end
    end
end

function err = cost_function_entropy_scaled(x, cases, n1, Q, R, T, Sc_ratio, Dfc)
    Sg = x(1); B_vec = x(2:end); exp_L = [cases.exp_life];
    penalty = 0; p_lives = zeros(1,5);
    
    if Sg <= 0 || any(B_vec <= 0), err = 1e12; return; end
    
    for k = 1:5
        % CORRECTED: Match the 4 outputs of the helper function
        [Dc, Delta_Sf, ~, ~] = calculate_advanced_metrics(cases(k), B_vec(k), n1, Q, R, T, k);
        
        den = log(1 - Sc_ratio); 
        Df = max(1e-12, ((Dfc / den) * log(max(1e-12, 1 - Delta_Sf / Sg))));
        
        % Force Df to be 1 order of magnitude larger than Dc
        % If Df/Dc is not around 10, add to penalty
        target_ratio = 10;
        current_ratio = Df / max(1e-15, Dc);
        penalty = penalty + (log10(current_ratio) - log10(target_ratio))^2 * 50;
        
        p_lives(k) = 1 / (Dc + Df);
    end
    err = sqrt(mean((log10(exp_L) - log10(p_lives)).^2)) + penalty;
end