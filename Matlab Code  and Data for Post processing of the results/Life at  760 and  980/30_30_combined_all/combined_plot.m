clc;
clear;
close all;

%% 1. SETUP: DEFINE THE FOUR CASES (30/30 Combined Dwell)
% =========================================================================
cases(1).filename = 'shi_30_30_2_4.csv'; 
cases(1).epsilon_a = 0.012;      
cases(1).exp_life = 5.6;           

cases(2).filename = 'pin_30_30_2_0.csv'; 
cases(2).epsilon_a = 0.010;      
cases(2).exp_life = 170;         

cases(3).filename = 'Fatigue.csv'; 
cases(3).epsilon_a = 0.009;      
cases(3).exp_life = 370;          

cases(4).filename = 'shi_30_30_1_6.csv'; 
cases(4).epsilon_a = 0.008;      
cases(4).exp_life = 562;         

%% GLOBAL MATERIAL CONSTANTS
Sg = 0.38;      
Sc = 0.9*Sg;   
Dfc = 0.9; 
n1 = 0.6;        
Q = 6.97e-19;    
T = 1033;        
R = 8.314;       
q_int = 0.4; % Interaction exponent

% Pre-allocate storage for LDS analysis
pred_lives_linear = zeros(1,4);
pred_lives_q      = zeros(1,4);
exp_lives         = zeros(1,4);
Dc_values         = zeros(1,4); % NEW: Creep Damage Storage
Df_values         = zeros(1,4); % NEW: Fatigue Damage Storage

%% 2. MAIN LOOP
for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    
    % B1 coefficient based on strain amplitude
  B1 = 2.0428e3 + (-2.1056e4 * eps_a) + (-1.2398e7 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:,1);
        entropy_rate_col = data(:,25);
        
        Delta_Sf = 0; 
        Dc = 0;
        
        for i = 2:length(time_col)
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            t_curr = time_col(i);
            
            is_creep = false; 
            is_fatigue = false;
            
            % ----- CASE TIME WINDOWS (30s/30s Dwell Logic) -----
            if k == 1
                if (t_curr>660 && t_curr<=690)||(t_curr>714 && t_curr<=744), is_creep = true;
                elseif (t_curr>690 && t_curr<=714)||(t_curr>744 && t_curr<=768), is_fatigue = true; end
            elseif k == 2
                if (t_curr>10 && t_curr<=40)||(t_curr>60 && t_curr<=90), is_creep = true;
                elseif (t_curr>40 && t_curr<=60)||(t_curr>90 && t_curr<=110), is_fatigue = true; end
            elseif k == 3
                if (t_curr>681 && t_curr<=711)||(t_curr>729 && t_curr<=759), is_creep = true;
                elseif (t_curr>711 && t_curr<=729)||(t_curr>759 && t_curr<=777), is_fatigue = true; end
            elseif k == 4
                if (t_curr>560 && t_curr<=590)||(t_curr>606 && t_curr<=636), is_creep = true;
                elseif (t_curr>590 && t_curr<=606)||(t_curr>636 && t_curr<=652), is_fatigue = true; end
            end
            
            % Accumulate Damage
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + s_dot*dt;
            end
        end
        
        % Calculate Final Fatigue Damage Df
        if Delta_Sf >= Sg
            Df = 1;
        else
            Df = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
        end
        
        % Store LDS Damage Values
        Dc_values(k) = Dc;
        Df_values(k) = Df;
        
        % Linear Prediction
        pred_lives_linear(k) = 1/(Dc + Df);
        % Interaction Prediction (q=0.4)
        pred_lives_q(k) = 1/(Dc^q_int + Df^q_int);
    else
        pred_lives_linear(k) = NaN; 
        pred_lives_q(k) = NaN;
        Dc_values(k) = NaN;
        Df_values(k) = NaN;
    end
    exp_lives(k)  = cases(k).exp_life;
end

%% 3. OUTPUT DAMAGE SUMMARY TABLE
% =========================================================================
fprintf('\n--- Damage Summation Results (LDS) ---\n');
LDS_Table = table([1:4]', [cases.epsilon_a]', exp_lives', Dc_values', Df_values', ...
    'VariableNames', {'Case', 'Strain_Amp', 'Exp_Life', 'Dc_per_cycle', 'Df_per_cycle'});
disp(LDS_Table);

%% 4. CORRELATION PLOTS
% =========================================================================
createCorrelationPlot(exp_lives, pred_lives_linear, 'Linear Summation', 'N_f = (D_c + D_f)^{-1}');
createCorrelationPlot(exp_lives, pred_lives_q, ['Interaction (q=', num2str(q_int), ')'], ...
    ['N_f = (D_c^{', num2str(q_int), '} + D_f^{', num2str(q_int), '})^{-1}']);

%% 5. STRAIN RANGE vs LIFE COMPARISON
% =========================================================================
strain_range = 2 * [cases.epsilon_a];
figure('Color','w','Position', [100 100 800 600]); hold on;
hExp = plot(exp_lives, strain_range, 'ok', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
hLin = plot(pred_lives_linear, strain_range, 'sr', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
hInt = plot(pred_lives_q, strain_range, 'db', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

% Smooth trendline for the interaction model
valid_idx = ~isnan(pred_lives_q);
p = polyfit(log10(pred_lives_q(valid_idx)), strain_range(valid_idx), 2);
x_smooth = logspace(log10(min(pred_lives_q(valid_idx))), log10(max(pred_lives_q(valid_idx))), 100);
y_smooth = polyval(p, log10(x_smooth));
hFit = plot(x_smooth, y_smooth, 'b--', 'LineWidth', 2);

set(gca, 'XScale', 'log', 'FontSize', 18, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlabel('Cycles to Failure, N_f', 'FontSize', 22);
ylabel('Strain Range, \Delta\epsilon', 'FontSize', 22);
grid on; box on;
legend([hExp, hLin, hInt, hFit], {'Experiment', 'Linear Sum', 'Interaction Model', 'Sim. Trendline'}, 'Location', 'northeast');
title('Strain-Life Comparison: DD6 at 760°C');
hold off;

%% HELPER FUNCTION
function createCorrelationPlot(exp, pred, pTitle, formula)
    figure('Color','w'); hold on;
    colors = {'r','b','g','m'}; markers = {'s','d','^','v'};
    all_v = [exp, pred];
    lims = [min(all_v)*0.1, max(all_v)*10];
    
    loglog(lims, lims, 'k-', 'LineWidth', 2);
    loglog(lims, lims*2, 'k--', 'LineWidth', 1.5); loglog(lims, lims/2, 'k--', 'LineWidth', 1.5);
    loglog(lims, lims*3, 'b-.', 'LineWidth', 1.5); loglog(lims, lims/3, 'b-.', 'LineWidth', 1.5);
    
    for k = 1:length(exp)
        if ~isnan(pred(k))
            loglog(exp(k), pred(k), markers{k}, 'MarkerSize', 12, 'LineWidth', 1.5, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
        end
    end
    
    grid on; axis square; box on; xlim(lims); ylim(lims);
    set(gca, 'LineWidth', 1.5, 'FontSize', 18, 'FontName', 'Times New Roman');
    xlabel('Experimental Life (Cycles)'); ylabel('Predicted Life (Cycles)');
    title(pTitle);
    annotation('textbox', [0.55 0.15 0.35 0.1], 'String', formula, 'Interpreter', 'tex', ...
        'FontSize', 14, 'FontName', 'Times New Roman', 'BackgroundColor', 'white');
end