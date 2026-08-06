clc;
clear;
close all;

%% ========================================================================
%% 1. SETUP: DEFINE THE CASES
%% ========================================================================

cases(1).filename = '0_0_R__1.csv'; 
cases(1).epsilon_a = 0.01;      
cases(1).exp_life = 5;           

cases(2).filename = 'mean_strain_0_002.csv'; 
cases(2).epsilon_a = 0.01;      
cases(2).exp_life = 600;         

cases(3).filename = 'mean_strain_0_01.csv'; 
cases(3).epsilon_a = 0.01;      
cases(3).exp_life = 297;          

cases(4).filename = 'mean_strain_0_015.csv'; 
cases(4).epsilon_a = 0.01;      
cases(4).exp_life = 100;         

numCases = length(cases);

%% ========================================================================
%% 2. GLOBAL MATERIAL CONSTANTS
%% ========================================================================

n1 = 0.6;        
Q  = 6.97e-19;    
T  = 1033;        
R  = 8.314;       
Sg = 0.38;      
Sc = 0.9*Sg;   
Dfc = 0.9;  
q = 0.576;    

% Initialize arrays
pred_lives_linear = zeros(1,numCases);
pred_lives_q      = zeros(1,numCases);
exp_lives         = zeros(1,numCases);

%% ========================================================================
%% 3. MAIN COMPUTATION LOOP
%% ========================================================================

for k = 1:numCases
    
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        
        data = readmatrix(current_file);
        time_col = data(:,1);
        entropy_rate_col = data(:,25);
        
        Delta_Sf = 0; 
        Dc = 0;
        
        for i = 2:length(time_col)
            
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false; 
            is_fatigue = false;
            
            % ----- CASE TIME WINDOWS -----
            if k==1
                if (t_current>0 && t_current<=0), is_creep=true;
                elseif (t_current>10 && t_current<=50), is_fatigue=true; end
                
            elseif k==2
                if (t_current>12 && t_current<=72), is_creep=true;
                elseif (t_current>72 && t_current<=112), is_fatigue=true; end
                
            elseif k==3
                if (t_current>20 && t_current<=80), is_creep=true;
                elseif (t_current>80&& t_current<=120), is_fatigue=true; end
                
            elseif k==4
                if (t_current>25 && t_current<=85), is_creep=true;
                elseif (t_current>85 && t_current<=125), is_fatigue=true; end
            end
            
            % ---- Creep Damage ----
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi) * (s_dot)^(1-n1) * dt;
            end
            
            % ---- Fatigue Damage ----
            if is_fatigue
                Delta_Sf = Delta_Sf + s_dot * dt;
            end
        end
        
        % ---- Fatigue Damage Model ----
        if Delta_Sf >= Sg
            Df = 1;
        else
            Df = -Dfc/log(1-Sc/Sg) * log(1-Delta_Sf/Sg);
        end
        
        % ---- Life Predictions ----
        pred_lives_linear(k) = 1 / (Dc + Df);
        pred_lives_q(k)      = 1 / (Dc^q + Df^q);
        
    else
        pred_lives_linear(k) = NaN;
        pred_lives_q(k)      = NaN;
    end
    
    exp_lives(k)  = cases(k).exp_life;
end

%% ========================================================================
%% 4. PLOT 1: LINEAR SUMMATION MODEL
%% ========================================================================

plotCorrelation(exp_lives, pred_lives_linear, ...
    'Linear Summation Model', ...
    'N_f = (D_c + D_f)^{-1}');

%% ========================================================================
%% 5. PLOT 2: INTERACTION MODEL
%% ========================================================================

plotCorrelation(exp_lives, pred_lives_q, ...
    ['Interaction Model (q = ', num2str(q), ')'], ...
    ['N_f = (D_c^{', num2str(q), '} + D_f^{', num2str(q), '})^{-1}']);

%% ========================================================================
%% 6. PLOT 3: STRAIN RANGE vs LIFE
%% ========================================================================

strain_range = 2 * [cases.epsilon_a];

figure('Color','w','Position', [100 100 800 800]); 
hold on;

hExp = plot(exp_lives, strain_range, 'ok', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');

hLin = plot(pred_lives_linear, strain_range, 'sr', ...
    'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

hInt = plot(pred_lives_q, strain_range, 'db', ...
    'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

set(gca, 'XScale', 'log', ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'FontName', 'Times New Roman');

xlabel('Cycles to Failure, N_f', 'FontSize', 26);
ylabel('Strain Range, \Delta\epsilon', 'FontSize', 26);

grid on; 
box on; 
axis square;

legend([hExp, hLin, hInt], ...
    {'Experiment', 'Linear Model', 'Interaction Model'}, ...
    'Location', 'northeast');

title('Model Comparison', 'FontWeight', 'normal');

hold off;

%% ========================================================================
%% HELPER FUNCTION
%% ========================================================================

function plotCorrelation(x, y, plotTitle, formulaStr)

    figure('Color','w'); 
    hold on;
    
    colors  = lines(length(x)); 
    markers = {'s','d','^','v','p','o','>'};
    
    all_vals = [x(:); y(:)];
    limits = [min(all_vals)*0.5, max(all_vals)*5];
    
    % Perfect and tolerance lines
    loglog(limits, limits, 'k-', 'LineWidth', 2);
    loglog(limits, limits*2, 'k--', 'LineWidth', 1.5);
    loglog(limits, limits*0.5, 'k--', 'LineWidth', 1.5);
    
    % Data points
    for k = 1:length(x)
        loglog(x(k), y(k), markers{k}, ...
            'MarkerSize', 12, ...
            'LineWidth', 1.5, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(k,:));
    end
    
    grid on; 
    axis square; 
    box on;
    
    xlim(limits); 
    ylim(limits);
    
    set(gca, 'LineWidth', 1.5, ...
        'FontSize', 22, ...
        'FontName', 'Times New Roman');
    
    xlabel('Experimental Life (Cycles)');
    ylabel('Predicted Life (Cycles)');
    title(plotTitle, 'FontWeight', 'normal');
    
    annotation('textbox', [0.6 0.15 0.3 0.1], ...
        'String', {formulaStr}, ...
        'Interpreter', 'tex', ...
        'FontSize', 18, ...
        'FontName', 'Times New Roman', ...
        'BackgroundColor', 'white');
end
