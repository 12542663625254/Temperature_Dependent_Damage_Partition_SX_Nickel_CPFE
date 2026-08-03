clc;
clear;
close all;

%% 1. SETUP: DEFINE THE FIVE CASES (0/60 Dwell)
% =========================================================================
% CASE 1
cases(1).filename = 'Pin_2_4_60C.csv'; 
cases(1).epsilon_a = 0.012;      
cases(1).exp_life = 16;           
% CASE 2
cases(2).filename = 'Pin_2_0_60C.csv'; 
cases(2).epsilon_a = 0.01;      
cases(2).exp_life = 129;         
% CASE 3
cases(3).filename = 'Pin_1_8_60C.csv'; 
cases(3).epsilon_a = 0.009;      
cases(3).exp_life = 111;          
% CASE 4
cases(4).filename = 'Pin_1_6_60_Comp.csv'; 
cases(4).epsilon_a = 0.008;      
cases(4).exp_life = 473;         
% CASE 5
cases(5).filename = '0_60_0_007.csv';   
cases(5).epsilon_a = 0.007;             
cases(5).exp_life = 1000;               

% GLOBAL CONSTANTS
n1 = 0.6;        
Q = 6.97e-19;    
T = 1033;        
R = 8.314;       
D0 = 0; 

%% 2. PRE-CALCULATE CREEP DAMAGE (Dc) & DELTA_SF
% =========================================================================
fprintf('Reading files and pre-calculating creep/fatigue inputs...\n');

Dc_stored = zeros(1, 5);
Delta_Sf_stored = zeros(1, 5);

for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    
    % B1 and PHI Calculation
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:,1);
        entropy_rate_col = data(:,25);
        
        Delta_Sf_accum = 0; 
        Dc_accum = 0;       
        num_points = length(time_col);
        
        for i = 2:num_points
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false; is_fatigue = false;
            
            % --- 0/60 TIMING LOGIC ---
            if k == 1
                if (t_current > 252 && t_current <= 312)
                    is_creep = true;
                elseif (t_current > 312 && t_current <= 336) || (t_current > 228 && t_current <= 252)
                    is_fatigue = true;
                end
            elseif k == 2
                if (t_current > 330 && t_current <= 390)
                    is_creep = true;
                elseif (t_current > 390 && t_current <= 410) || (t_current > 310 && t_current <= 330)
                    is_fatigue = true;
                end
            elseif k == 3
                if (t_current > 315 && t_current <= 375)
                    is_creep = true;
                elseif (t_current > 297 && t_current <= 315) || (t_current > 375 && t_current <= 393)
                    is_fatigue = true;
                end
            elseif k == 4
                if (t_current > 116 && t_current <= 176)
                    is_creep = true;
                elseif (t_current > 176 && t_current <= 192) || (t_current > 100 && t_current <= 116)
                    is_fatigue = true;
                end
            elseif k == 5
                if (t_current > 81 && t_current <= 141) 
                    is_creep = true;
                elseif (t_current > 141 && t_current <= 155) || (t_current > 67 && t_current <= 81)
                    is_fatigue = true;
                end
            end
            
            if is_creep && s_dot > 0
                Dc_accum = Dc_accum + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf_accum = Delta_Sf_accum + s_dot*dt;
            end
        end
        Dc_stored(k) = Dc_accum;
        Delta_Sf_stored(k) = Delta_Sf_accum;
    else
        Dc_stored(k) = NaN; Delta_Sf_stored(k) = NaN;
    end
end

%% 3. RUN OPTIMIZATION (fminsearch)
% =========================================================================
% Initial Guess [Sg, Sc_ratio, Dfc]
x0 = [0.013, 0.95, 0.1]; 

options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 2000);
exp_lives = [cases.exp_life];
objFun = @(x) cost_function(x, Dc_stored, Delta_Sf_stored, D0, exp_lives);

fprintf('\nStarting Optimization for 0/60 Dwell...\n');
[x_opt, fval] = fminsearch(objFun, x0, options);

Sg_opt = x_opt(1);
Sc_ratio_opt = x_opt(2);
Dfc_opt = x_opt(3);
Sc_opt = Sc_ratio_opt * Sg_opt;

fprintf('\n==============================================\n');
fprintf('OPTIMIZATION RESULTS (0/60 Data):\n');
fprintf('==============================================\n');
fprintf('Sg (Optimized)   = %.5f\n', Sg_opt);
fprintf('Sc (Optimized)   = %.5f\n', Sc_opt);
fprintf('Sc/Sg Ratio      = %.5f\n', Sc_ratio_opt);
fprintf('Dfc (Optimized)  = %.5f\n', Dfc_opt);
fprintf('Final Residual   = %.5f\n', fval);
fprintf('==============================================\n');

% COPY THESE VALUES INTO YOUR PLOTTING SCRIPT!

%% ========================================================================
%% COST FUNCTION
%% ========================================================================
function err = cost_function(x, Dc_vec, Delta_Sf_vec, D0, exp_lives)
    Sg_curr = x(1);
    Sc_ratio = x(2);
    Dfc_curr = x(3);
    
    penalty = 0;
    if Sg_curr <= 0.001 || Sg_curr > 0.5, penalty = penalty + 1e5; end
    if Sc_ratio <= 0.01 || Sc_ratio >= 0.999, penalty = penalty + 1e5; end
    if Dfc_curr <= 0.001 || Dfc_curr > 10, penalty = penalty + 1e5; end
    
    if penalty > 0, err = 1e6 + penalty; return; end

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
            
            if val_inside <= 1e-9
                term_log = -50; 
            else
                term_log = log(val_inside);
            end
            
            Df = D0 + (numerator / denominator) * term_log;
        end
        if Df < 0, Df = 0; end
        
        total = Dc + Df;
        if total <= 1e-9, pred_lives(k) = 1e9; else, pred_lives(k) = 1/total; end
    end
    
    log_exp = log10(exp_lives);
    log_pred = log10(pred_lives);
    err = sqrt(mean((log_exp - log_pred).^2)) + penalty;
end