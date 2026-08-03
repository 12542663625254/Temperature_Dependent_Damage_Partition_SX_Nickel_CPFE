clc; clear; close all;

%% 1. SETUP: DEFINE THE FIVE CASES (DD6 Superalloy @ 760°C)
% =========================================================================
cases(1).filename = 'shi_30_30_2_4.csv'; cases(1).epsilon_a = 0.012; cases(1).exp_life = 5;           
cases(2).filename = 'pin_30_30_2_0.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 250;         
cases(3).filename = 'Fatigue.csv';       cases(3).epsilon_a = 0.009; cases(3).exp_life = 370;          
cases(4).filename = 'shi_30_30_1_6.csv'; cases(4).epsilon_a = 0.008; cases(4).exp_life = 562;         
cases(5).filename = '30_30_0_007.csv';   cases(5).epsilon_a = 0.007; cases(5).exp_life = 1455;                 

% GLOBAL CONSTANTS
n1 = 0.6; Q = 6.97e-19; T = 1033; R = 8.314; D0 = 0; 

%% 2. FIXED MODEL PARAMETERS
% =========================================================================
Sc_ratio_fixed = 0.9;
Dfc_fixed = 0.9;

%% 3. DATA PRE-LOADING & OPTIMIZATION
% =========================================================================
% x0: [Initial Sg, B1_1, B1_2, B1_3, B1_4, B1_5]
x0 = [0.8, 3.8e5, 3.8e5, 3.8e5, 3.8e5, 3.8e5]; 
options = optimset('Display', 'iter', 'MaxFunEvals', 40000, 'MaxIter', 40000);
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

objFun = @(x) cost_function_Sg_B_variable(x, cases, n1, Q, R, T, D0, ...
              exp_lives, Sc_ratio_fixed, Dfc_fixed);

fprintf('\nStarting Optimization for Sg and Strain-Specific B1 values...\n');
[x_opt, fval] = fminsearch(objFun, x0, options);

Sg_opt = x_opt(1);
B1_opt_vec = x_opt(2:end);

%% 4. FINAL DAMAGE CALCULATION (FOR TABLE)
% =========================================================================
pred_lives = zeros(1,5);
Dc_vals = zeros(1,5);
Df_vals = zeros(1,5);

for k = 1:5
    [Dc, Delta_Sf] = calculate_case_damage(cases(k), B1_opt_vec(k), n1, Q, R, T, k);
    
    if Delta_Sf >= Sg_opt
        Df = 1.0;
    else
        den = log(1 - Sc_ratio_fixed); 
        val_in = 1 - Delta_Sf / Sg_opt;
        log_t = log(max(1e-10, val_in));
        Df = D0 + ((Dfc_fixed - D0) / den) * log_t;
    end
    
    Dc_vals(k) = Dc;
    Df_vals(k) = max(0, Df);
    pred_lives(k) = 1 / (Dc_vals(k) + Df_vals(k));
end

%% 5. POLYNOMIAL FITTING (B1 vs Strain Amplitude)
% =========================================================================
strain_amps = [cases.epsilon_a]';
B1_vals = B1_opt_vec';

% Fit a 2nd-degree polynomial: B1 = p1*ea^2 + p2*ea + p3
p_coeffs = polyfit(strain_amps, B1_vals, 2);
B1_poly_func = @(ea) polyval(p_coeffs, ea);

ea_range = linspace(min(strain_amps)*0.9, max(strain_amps)*1.1, 100);
B1_poly_curve = B1_poly_func(ea_range);

%% 6. DISPLAY RESULTS
% =========================================================================
fprintf('\n============================================================\n');
fprintf('     FINAL OPTIMIZED RESULTS (Sg and B_1 varied)           \n');
fprintf('============================================================\n');
fprintf('Optimized Fatigue Limit (Sg):      %.6f\n', Sg_opt);
fprintf('Polynomial Coeffs (B1=p1*ea^2+p2*ea+p3):\n');
fprintf('  p1: %.4e, p2: %.4e, p3: %.4e\n', p_coeffs(1), p_coeffs(2), p_coeffs(3));
fprintf('------------------------------------------------------------\n');

Error_Pct = abs(pred_lives - exp_lives) ./ exp_lives * 100;
T_results = table((1:5)', strain_amps, B1_vals, exp_lives', pred_lives', ...
    Dc_vals', Df_vals', Error_Pct', ...
    'VariableNames', {'Case', 'Strain_Amp', 'Opt_B1', 'Exp_Life', 'Pred_Life', 'Dc', 'Df', 'Err_Pct'});
disp(T_results);

%% 7. VISUALIZATION (Including Damage Interaction Plot)
% =========================================================================
colors = {'r', 'b', 'g', 'm', 'k'};

% --- FIGURE 1: B1 Polynomial Fit ---
figure('Color', 'w', 'Position', [50, 100, 500, 450]);
plot(strain_amps, B1_vals, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Optimized B_1');
hold on;
plot(ea_range, B1_poly_curve, 'm-', 'LineWidth', 2, 'DisplayName', '2nd-Order Poly Fit');
grid on; set(gca, 'FontSize', 11, 'TickDir', 'in');
xlabel('Strain Amplitude (\epsilon_a)'); ylabel('Creep Parameter (B_1)');
title('B_1 Correlation'); legend('Location', 'best');

% --- FIGURE 2: Prediction Accuracy ---
figure('Color', 'w', 'Position', [560, 100, 500, 450]);
loglog([min(exp_lives)/2, max(exp_lives)*2], [min(exp_lives)/2, max(exp_lives)*2], 'k-', 'LineWidth', 2); hold on;
loglog([min(exp_lives)/2, max(exp_lives)*2], [min(exp_lives)/2, max(exp_lives)*2]*2, 'k--');
loglog([min(exp_lives)/2, max(exp_lives)*2], [min(exp_lives)/2, max(exp_lives)*2]*0.5, 'k--');
for k = 1:5
    loglog(exp_lives(k), pred_lives(k), 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors{k}, 'DisplayName', sprintf('ea=%.3f', strain_amps(k)));
end
grid on; axis square; xlabel('Exp Life (Cycles)'); ylabel('Pred Life (Cycles)'); title('Life Prediction Accuracy');

% --- FIGURE 3: CREEP-FATIGUE DAMAGE INTERACTION ---
figure('Color', 'w', 'Position', [1070, 100, 550, 450]);
% Plot the Linear Damage Rule Line (Dc + Df = 1)
df_line = linspace(0, 1, 100);
dc_line = 1 - df_line;
plot(df_line, dc_line, 'r-', 'LineWidth', 2, 'DisplayName', 'Equivalent Damage Line'); 
hold on;

% Plot the optimized damage points for each case
for k = 1:5
    % For a cycle-based interaction, we multiply damage by Exp_Life
    % to see where failure occurred on the Dc-Df map
    Df_total = Df_vals(k) * exp_lives(k);
    Dc_total = Dc_vals(k) * exp_lives(k);
    
    plot(Df_total, Dc_total, 'Marker', 's', 'MarkerSize', 12, ...
        'MarkerFaceColor', colors{k}, 'MarkerEdgeColor', 'k', ...
        'LineStyle', 'none', 'DisplayName', sprintf('Case %d', k));
end

grid on; axis square;
xlabel('Fatigue Damage (D_f)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Creep Damage (D_c)', 'FontSize', 12, 'FontWeight', 'bold');
title('Damage Interaction Diagram (DD6 @ 760°C)', 'FontSize', 13);
xlim([0 1.2]); ylim([0 1.2]);
legend('Location', 'northeast', 'FontSize', 9);
set(gca, 'TickDir', 'in', 'LineWidth', 1.2);
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
            if (t_curr>228 && t_curr<=258)||(t_curr>282&& t_curr<=312), is_creep = true;
            elseif (t_curr>258 && t_curr<=282)||(t_curr>312 && t_curr<=336), is_fatigue = true; end
        elseif k == 2
            if (t_curr>10 && t_curr<=40)||(t_curr>60 && t_curr<=90), is_creep = true;
            elseif (t_curr>40 && t_curr<=60)||(t_curr>90 && t_curr<=110), is_fatigue = true; end
        elseif k == 3
            if (t_curr>9 && t_curr<=39)||(t_curr>57 && t_curr<=87), is_creep = true;
            elseif (t_curr>39 && t_curr<=57)||(t_curr>87 && t_curr<=105), is_fatigue = true; end
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

function err = cost_function_Sg_B_variable(x, cases, n1, Q, R, T, D0, exp_lives, Sc_ratio, Dfc)
    Sg_curr = x(1); B_vec = x(2:end); p_lives = zeros(1,5);
    if Sg_curr <= 0 || any(B_vec <= 0), err = 1e10; return; end
    for k = 1:5
        [Dc, Delta_Sf] = calculate_case_damage(cases(k), B_vec(k), n1, Q, R, T, k);
        if Delta_Sf >= Sg_curr
            Df = 1.0;
        else
            den = log(1 - Sc_ratio); val_in = 1 - Delta_Sf / Sg_curr;
            log_t = log(max(1e-10, val_in));
            Df = D0 + ((Dfc - D0) / den) * log_t;
        end
        p_lives(k) = 1 / (Dc + max(0, Df));
    end
    err = sqrt(mean((log10(exp_lives) - log10(p_lives)).^2));
end