clc; clear; close all;

%% 1. SETUP: DATA & CONSTANTS
% =========================================================================
cases(1).filename = 'shi_30_30_2_4.csv'; cases(1).epsilon_a = 0.012; cases(1).exp_life = 30;           
cases(2).filename = 'shi_30_30_2_0.csv'; cases(2).epsilon_a = 0.010; cases(2).exp_life = 290;         
cases(3).filename = 'shi_30_30_1_8.csv'; cases(3).epsilon_a = 0.009; cases(3).exp_life = 850;          
cases(4).filename = 'shi_30_30_1_6.csv'; cases(4).epsilon_a = 0.008; cases(4).exp_life = 3557;         

% GLOBAL CONSTANTS (FIXED)
n1 = 0.6; Q = 6.97e-19; T = 1033; R = 8.314; D0 = 0; 

%% 2. PRE-CALCULATE DAMAGE EVOLUTION (CREEP-FATIGUE INTERACTION)
% =========================================================================
fprintf('Reading files and tracking damage interaction...\n');
Dc_stored = zeros(1, 4);
Delta_Sf_stored = zeros(1, 4);

% Storage for evolution plotting
time_history = cell(1, 4);
Dc_history = cell(1, 4);
Sf_history = cell(1, 4);

for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        t_vec = data(:, 1);          
        s_dot_vec = data(:, 25); 
        
        Dc_evol = zeros(size(t_vec));
        Sf_evol = zeros(size(t_vec));
        
        for i = 2:length(t_vec)
            dt = t_vec(i) - t_vec(i-1);
            s_dot = s_dot_vec(i);
            t_curr = t_vec(i);
            
            is_fatigue = false;
            is_creep = false;
            
            % --- UPDATED TIME LOGIC ---
            if k == 1
                if (t_curr > 80 && t_curr <= 120),    is_fatigue = true;
                elseif (t_curr > 65 && t_curr <= 80), is_creep = true; end
            elseif k == 2
                if (t_curr > 80 && t_curr <= 120),    is_fatigue = true;
                elseif (t_curr > 65 && t_curr <= 80), is_creep = true; end
            elseif k == 3
                if (t_curr > 80 && t_curr <= 120),    is_fatigue = true;
                elseif (t_curr > 65 && t_curr <= 80), is_creep = true; end
            elseif k == 4
                if (t_curr > 80 && t_curr <= 120),    is_fatigue = true;
                elseif (t_curr > 65 && t_curr <= 80), is_creep = true; end
            end
            
            % Accumulate Creep Damage
            if is_creep && s_dot > 0
                term = (1/phi) * (s_dot)^(1 - n1);
                Dc_evol(i) = Dc_evol(i-1) + (term * dt);
            else
                Dc_evol(i) = Dc_evol(i-1);
            end
            
            % Accumulate Fatigue Entropy (Delta_Sf)
            if is_fatigue
                Sf_evol(i) = Sf_evol(i-1) + (s_dot * dt);
            else
                Sf_evol(i) = Sf_evol(i-1);
            end
        end
        
        % Store histories and final cycle values
        time_history{k} = t_vec;
        Dc_history{k} = Dc_evol;
        Sf_history{k} = Sf_evol;
        Dc_stored(k) = Dc_evol(end);
        Delta_Sf_stored(k) = Sf_evol(end);
    end
end

%% 3. OPTIMIZATION SETUP
% =========================================================================
x0 = [0.038, 0.9, 0.9]; 
options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 2000);
exp_lives = [cases.exp_life];
objFun = @(x) cost_function(x, Dc_stored, Delta_Sf_stored, D0, exp_lives);

[x_opt, fval] = fminsearch(objFun, x0, options);

Sg_opt = x_opt(1); Sc_opt = x_opt(2)*x_opt(1); Dfc_opt = x_opt(3);

%% 4. RESULTS & PLOTTING
% =========================================================================
pred_lives = zeros(1,4);
for k = 1:4
    Df_final = calc_Df(Delta_Sf_stored(k), Sg_opt, Sc_opt, Dfc_opt, D0);
    total_damage = Dc_stored(k) + Df_final;
    pred_lives(k) = 1 / total_damage;
end

% PLOT 1: Accuracy
figure('Color', 'w', 'Name', 'Life Prediction');
loglog([10 100000], [10 100000], 'k-', 'LineWidth', 2); hold on;
loglog([10 100000], [20 200000], 'k--', 'LineWidth', 1);
loglog([10 100000], [5 50000], 'k--', 'LineWidth', 1);
colors = {'r', 'b', 'g', 'm'}; markers = {'o', 's', '^', 'd'};
for k = 1:4
    loglog(exp_lives(k), pred_lives(k), markers{k}, 'MarkerSize', 12, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
end
grid on; axis square; xlabel('Experimental Life (N_f)'); ylabel('Predicted Life (N_f)');
title('Creep-Fatigue Interaction Fitting');

% PLOT 2: Damage Evolution

figure('Color', 'w', 'Name', 'Damage History', 'Position', [100 100 1100 800]);
for k = 1:4
    subplot(2, 2, k); hold on;
    Df_evol = arrayfun(@(sf) calc_Df(sf, Sg_opt, Sc_opt, Dfc_opt, D0), Sf_history{k});
    
    plot(time_history{k}, Dc_history{k}, 'r--', 'LineWidth', 2, 'DisplayName', 'D_c (Creep)');
    plot(time_history{k}, Df_evol, 'b-', 'LineWidth', 2, 'DisplayName', 'D_f (Fatigue)');
    plot(time_history{k}, Dc_history{k} + Df_evol, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Total');
    
    title(sprintf('Case %d (\\epsilon_a=%.3f)', k, cases(k).epsilon_a));
    xlabel('Time (s)'); ylabel('Damage');
    if k==1, legend('Location', 'northwest', 'FontSize', 8); end
    grid on; box on;
end

%% HELPER FUNCTIONS
function Df = calc_Df(Sf, Sg, Sc, Dfc, D0)
    if Sf >= Sg, Df = 1.0;
    else
        den = log(1 - Sc/Sg);
        val = 1 - Sf/Sg;
        log_term = (val <= 1e-6) ? -20 : log(val);
        Df = D0 + ((Dfc - D0) / den) * log_term;
    end
    Df = max(0, Df);
end

function err = cost_function(x, Dc_v, Sf_v, D0, exp_v)
    if x(1)<=0 || x(2)<=0.01 || x(2)>=0.99 || x(3)<=0.1, err=1e9; return; end
    Sg = x(1); Sc = x(2)*x(1); Dfc = x(3);
    p_lives = zeros(size(exp_v));
    for k = 1:length(exp_v)
        total_D = Dc_v(k) + calc_Df(Sf_v(k), Sg, Sc, Dfc, D0);
        p_lives(k) = (total_D <= 1e-9) ? 1e9 : 1 / total_D;
    end
    err = sqrt(mean((log10(exp_v) - log10(p_lives)).^2));
end