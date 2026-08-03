clc;
clear;
close all;

%% 1. SETUP: DATA & CONSTANTS
% =========================================================================
% CASE 1: 1.2% Strain 
cases(1).filename = 'cyclic_2_4.csv'; 
cases(1).epsilon_a = 0.012;      
cases(1).exp_life = 30;           

% CASE 2: 0.8% Strain 
cases(2).filename = 'cyclic_2_0.csv'; 
cases(2).epsilon_a = 0.010;      
cases(2).exp_life = 290;         

% CASE 3: 0.9% Strain 
cases(3).filename = 'cyclic_1_8.csv'; 
cases(3).epsilon_a = 0.009;      
cases(3).exp_life = 850;          

% CASE 4: 1.0% Strain 
cases(4).filename = 'cyclic_1_6.csv'; 
cases(4).epsilon_a = 0.008;      
cases(4).exp_life = 3557;         

% GLOBAL CONSTANTS (FIXED)
n1 = 0.6;        
Q = 6.97e-19;    
T = 1033;        
R = 8.314;       
D0 = 0; % Initial damage

%% 2. PRE-CALCULATE CREEP DAMAGE (Dc) & DELTA_SF INPUTS
% =========================================================================
fprintf('Reading files and pre-calculating creep/fatigue inputs...\n');

Dc_stored = zeros(1, 4);
Delta_Sf_stored = zeros(1, 4);

for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    
    % B1 and PHI Calculation
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:, 1);          
        entropy_rate_col = data(:, 25); 
        
        Delta_Sf_accum = 0; 
        Dc_accum = 0;       
        num_points = length(time_col);
        
        for i = 2:num_points
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false;
            is_fatigue = false;
            
            % --- TIME LOGIC (PURE FATIGUE) ---
            if k == 1 
                if (t_current > 12 && t_current <= 60), is_fatigue = true; end
            elseif k == 2
               if (t_current > 10 && t_current <= 50), is_fatigue = true; end
            elseif k == 3
               if (t_current > 9 && t_current <= 45), is_fatigue = true; end
            elseif k == 4
                if (t_current > 8 && t_current <= 40), is_fatigue = true; end
            end
            
            % --- ACCUMULATION (DO NOT TOUCH) ---
            if is_creep && s_dot > 0
                term = (1/phi) * (s_dot)^(1 - n1);
                Dc_accum = Dc_accum + (term * dt);
            elseif is_fatigue
                Delta_Sf_accum = Delta_Sf_accum + (s_dot * dt);
            end
        end
        
        Dc_stored(k) = Dc_accum;
        Delta_Sf_stored(k) = Delta_Sf_accum;
    else
        % Fallback if file missing (optional warning)
        Dc_stored(k) = NaN;
        Delta_Sf_stored(k) = NaN;
    end
end

%% 3. OPTIMIZATION SETUP (Using fminsearch)
% =========================================================================
% Parameters: x = [Sg, Sc_ratio, Dfc]
% Initial Guesses
x0 = [0.038, 0.9, 0.9]; 

% Standard Optimization options
options = optimset('Display', 'iter', 'MaxFunEvals', 1000, 'MaxIter', 1000);

% Objective Function
exp_lives = [cases.exp_life];
objFun = @(x) cost_function(x, Dc_stored, Delta_Sf_stored, D0, exp_lives);

fprintf('\nStarting Optimization of Df parameters (using fminsearch)...\n');
[x_opt, fval] = fminsearch(objFun, x0, options);

% Extract Optimized Values
Sg_opt = x_opt(1);
Sc_ratio_opt = x_opt(2);
Dfc_opt = x_opt(3);
Sc_opt = Sc_ratio_opt * Sg_opt;

fprintf('\n==============================================\n');
fprintf('OPTIMIZATION RESULTS:\n');
fprintf('==============================================\n');
fprintf('Sg (Optimized)   = %.5f  (Original: 0.038)\n', Sg_opt);
fprintf('Sc (Optimized)   = %.5f  (Derived)\n', Sc_opt);
fprintf('Sc/Sg Ratio      = %.5f  (Original: 0.9)\n', Sc_ratio_opt);
fprintf('Dfc (Optimized)  = %.5f  (Original: 0.9)\n', Dfc_opt);
fprintf('Final Residual   = %.5f\n', fval);
fprintf('==============================================\n');

%% 4. CALCULATE FINAL LIVES WITH OPTIMIZED PARAMETERS
% =========================================================================
pred_lives_opt = zeros(1,4);

for k = 1:4
    Dc = Dc_stored(k);       
    Delta_Sf = Delta_Sf_stored(k);
    
    if Delta_Sf >= Sg_opt
        Df = 1.0;
    else
        numerator = Dfc_opt - D0;
        denominator = log(1 - Sc_opt/Sg_opt);
        % Protect log(negative)
        if denominator == 0, denominator = -1e-9; end
        
        % Protect against Delta_Sf > Sg inside log
        val_inside = 1 - Delta_Sf/Sg_opt;
        if val_inside <= 0
             log_term = -100; % Large damage
        else
             log_term = log(val_inside);
        end
        
        Df = D0 + (numerator / denominator) * log_term;
    end
    
    % Physics check: Damage cannot be negative
    if Df < 0, Df = 0; end
    
    total_damage = Dc + Df;
    if total_damage <= 0
        pred_lives_opt(k) = 1e8; 
    else
        pred_lives_opt(k) = 1 / total_damage;
    end
end

%% 5. PLOT RESULTS
% =========================================================================
figure('Color', 'w');
loglog([1 100000], [1 100000], 'k-', 'LineWidth', 2); hold on;
loglog([1 100000], [2 200000], 'k--', 'LineWidth', 1);
loglog([1 100000], [0.5 50000], 'k--', 'LineWidth', 1);

colors = {'r', 'b', 'g', 'm'};
markers = {'s', 'd', '^', 'v'};

for k = 1:4
    loglog(exp_lives(k), pred_lives_opt(k), markers{k}, ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
end

grid on; axis square;
xlabel('Experimental Life (N_f)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Predicted Life (N_f)', 'FontSize', 14, 'FontWeight', 'bold');
title('Optimization Result (No Toolbox)', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Perfect Match', '+2x Scatter', '-2x Scatter'}, 'Location', 'best');
xlim([10 100000]); ylim([10 100000]);


%% ========================================================================
%% HELPER FUNCTION: COST FUNCTION WITH PENALTIES
%% ========================================================================
function err = cost_function(x, Dc_vec, Delta_Sf_vec, D0, exp_lives)
    Sg_curr = x(1);
    Sc_ratio = x(2);
    Dfc_curr = x(3);
    
    % --- PENALTY FOR VIOLATING BOUNDS ---
    % Sg must be > 0
    % Sc_ratio must be between 0.01 and 0.99
    % Dfc must be > 0
    penalty = 0;
    if Sg_curr <= 0.001 || Sg_curr > 0.5, penalty = penalty + 1e5; end
    if Sc_ratio <= 0.01 || Sc_ratio >= 0.99, penalty = penalty + 1e5; end
    if Dfc_curr <= 0.1 || Dfc_curr > 10, penalty = penalty + 1e5; end
    
    % If parameters are totally invalid, return high error immediately
    if penalty > 0
        err = 1e6 + penalty;
        return;
    end

    Sc_curr = Sc_ratio * Sg_curr;
    pred_lives = zeros(size(exp_lives));
    
    for k = 1:length(exp_lives)
        Dc = Dc_vec(k);
        Delta_Sf = Delta_Sf_vec(k);
        
        if Delta_Sf >= Sg_curr
            Df = 1.0; 
        else
            numerator = Dfc_curr - D0;
            denominator = log(1 - Sc_curr / Sg_curr);
            
            val_inside = 1 - Delta_Sf / Sg_curr;
            if val_inside <= 1e-6
                term_log = -20; % represents failure approaching
            else
                term_log = log(val_inside);
            end
            
            Df = D0 + (numerator / denominator) * term_log;
        end
        
        if Df < 0, Df = 0; end
        
        total_damage = Dc + Df;
        if total_damage <= 1e-9
            pred_lives(k) = 1e9; 
        else
            pred_lives(k) = 1 / total_damage;
        end
    end
    
    % RMS Error of Log10 Life
    log_exp = log10(exp_lives);
    log_pred = log10(pred_lives);
    err = sqrt(mean((log_exp - log_pred).^2)) + penalty;
end