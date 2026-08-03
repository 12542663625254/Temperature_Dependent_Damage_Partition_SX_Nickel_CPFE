clc; clear; close all;

clc; clear; close all;

%% 1. SETUP: DEFINE THE CASES
% --- Existing 30/30 Cases ---
cases(1).epsilon_a = 0.012; cases(1).type = '30/30'; cases(1).filename = 'CFI_2_4_30_30_new.csv'; 
cases(2).epsilon_a = 0.011; cases(2).type = '30/30'; cases(2).filename = 'CFI_2_2_30_30.csv';
cases(3).epsilon_a = 0.010; cases(3).type = '30/30'; cases(3).filename = 'pin_30_30_2_0.csv';         
cases(4).epsilon_a = 0.009; cases(4).type = '30/30'; cases(4).filename = 'Fatigue_30_30_1_8.csv';          
cases(5).epsilon_a = 0.008; cases(5).type = '30/30'; cases(5).filename = 'shi_30_30_1_6.csv';         

% --- 60/0 Cases ---
cases(6).epsilon_a = 0.012; cases(6).type = '60/0';  cases(6).filename = 'Fatigue_60_0_2_4.csv'; 
cases(7).epsilon_a = 0.011; cases(7).type = '60/0';  cases(7).filename = 'CFI_2_2_60_0.csv'; 
cases(8).epsilon_a = 0.010; cases(8).type = '60/0';  cases(8).filename = 'Fatigue_60_0_2_0.csv'; 
cases(9).epsilon_a = 0.009; cases(9).type = '60/0';  cases(9).filename = 'Fatigue_60_0_1_8.csv'; 
cases(10).epsilon_a = 0.008; cases(10).type = '60/0'; cases(10).filename = 'Fatigue_60_0_0_8_perc.csv'; 

% --- 0/60 Cases ---
cases(11).epsilon_a = 0.012; cases(11).type = '0/60';  cases(11).filename = 'Pin_2_4_60C.csv'; 
cases(12).epsilon_a = 0.011; cases(12).type = '0/60';  cases(12).filename = 'CFI_2_2_0_60.csv'; 
cases(13).epsilon_a = 0.010; cases(13).type = '0/60';  cases(13).filename = 'Pin_2_0_60C.csv'; 
cases(14).epsilon_a = 0.009; cases(14).type = '0/60';  cases(14).filename = 'Pin_1_8_60C.csv'; 
cases(15).epsilon_a = 0.008; cases(15).type = '0/60';  cases(15).filename = 'CFI_0_8_0_60_new.csv'; 

% GLOBAL CONSTANTS
n1 = 0.6; D0 = 0;
Sg_fixed = 0.38; Sc_ratio_fixed = 0.9; Dfc_fixed = 0.9;
Q = 6.97e-19; T = 1033; R_univ = 8.314; 

Df_vals = zeros(length(cases),1);
Dc_vals = zeros(length(cases),1); 

%% 2. CALCULATION
for k = 1:length(cases)
    eps_k = cases(k).epsilon_a;
    B1_local = 2.0428e3 + (-2.1056e4 * eps_k) + (-1.2398e7 * eps_k^2);
   % B1_local = 389501 + (- 630849 * eps_k) + (257797 * eps_k^2);

    phi_local = B1_local * exp(-Q / (R_univ * T));
    
    if isfile(cases(k).filename)
        raw = readmatrix(cases(k).filename);
        t = raw(:, 1); s_dot_vec = raw(:, 25); 
        Delta_Sf = 0; C_sum = 0;
        
        for i = 2:length(t)
            dt = t(i) - t(i-1); s_dot = s_dot_vec(i);
            [is_c, is_f] = get_time_logic_expanded(t(i), k, cases(k).type);
            if is_c && s_dot > 0
                C_sum = C_sum + (1/phi_local) * (s_dot^(1-n1) * dt);
            elseif is_f
                Delta_Sf = Delta_Sf + (s_dot * dt);
            end
        end
        val_in = max(1e-10, 1 - Delta_Sf / Sg_fixed);
        Df_vals(k) = D0 + ((Dfc_fixed - D0) / log(1 - Sc_ratio_fixed)) * log(val_in);
        Dc_vals(k) = C_sum;
    else
        fprintf('File missing: %s\n', cases(k).filename);
        Df_vals(k) = NaN; Dc_vals(k) = NaN;
    end
end

%% 3. VISUALIZATION
fig = figure('Color','w','Position', [100 100 1100 800]);
hold on; box on; axis square;

% Logic indices
idx30   = find(strcmp({cases.type}, '30/30'));
idx60_0 = find(strcmp({cases.type}, '60/0'));
idx0_60 = find(strcmp({cases.type}, '0/60'));

% Color Definitions
color30   = [1, 0, 0];       % Red
color60_0 = [0, 0, 1];       % Blue
color0_60 = [0, 0.5, 0];     % Dark Green

% --- LEFT AXIS: Creep Damage (dc) ---
yyaxis left
p1 = plot([cases(idx30).epsilon_a], Dc_vals(idx30), '-s', 'Color', color30, 'LineWidth', 3, 'MarkerSize', 12, 'MarkerFaceColor', color30);
p2 = plot([cases(idx60_0).epsilon_a], Dc_vals(idx60_0), '-p', 'Color', color60_0, 'LineWidth', 2.5, 'MarkerSize', 14, 'MarkerFaceColor', color60_0);
p3 = plot([cases(idx0_60).epsilon_a], Dc_vals(idx0_60), '-o', 'Color', color0_60, 'LineWidth', 2.5, 'MarkerSize', 14, 'MarkerFaceColor', color0_60);

ylabel('Creep damage per cycle, $D_c$', 'Interpreter', 'latex', 'FontSize', 26);
set(gca, 'YScale', 'log', 'YColor', 'k'); 
ylim([1e-4, 1]);

% --- RIGHT AXIS: Fatigue Damage (df) ---
yyaxis right
p4 = plot([cases(idx30).epsilon_a], Df_vals(idx30), '--s', 'Color', color30, 'LineWidth', 3, 'MarkerSize', 12, 'MarkerFaceColor', color30);
p5 = plot([cases(idx60_0).epsilon_a], Df_vals(idx60_0), '--p', 'Color', color60_0, 'LineWidth', 2.5, 'MarkerSize', 14, 'MarkerFaceColor', color60_0);
p6 = plot([cases(idx0_60).epsilon_a], Df_vals(idx0_60), '--o', 'Color', color0_60, 'LineWidth', 2.5, 'MarkerSize', 14, 'MarkerFaceColor', color0_60);

ylabel('Fatigue damage per cycle, $D_f$', 'Interpreter', 'latex', 'FontSize', 26);
set(gca, 'YScale', 'log', 'YColor', 'k'); 
ylim([1e-4, 1]);

xlim([0.0079 0.0121]);
xlabel('Strain amplitude, $\epsilon_a$', 'Interpreter', 'latex', 'FontSize', 26);
grid on; 

% --- MAIN LEGEND (Updated to 6x1 layout) ---
legend([p1, p4, p2, p5, p3, p6], ...
    {'$D_c$ (30/30)', '$D_f$ (30/30)', '$D_c$ (60/0)', '$D_f$ (60/0)', '$D_c$ (0/60)', '$D_f$ (0/60)'}, ...
    'Location','northwest','Interpreter','latex','FontSize',22);

% --- ANNOTATION (Southwest Region - No Box) ---
% Position: [x_start y_start width height]
dim = [0.155 0.18 0.2 0.08]; 
str = {'Solid line: Creep damage ($D_c$)', 'Dashed line: Fatigue damage ($D_f$)'};

annotation('textbox', dim, 'String', str, 'Interpreter', 'latex', ...
    'FontSize', 20, ...
    'EdgeColor', 'none', ...       % CRITICAL: Removes the box border
    'BackgroundColor', 'none', ... % Set to 'none' to make it transparent
    'HorizontalAlignment', 'center', ...
    'FitBoxToText', 'on');

set(gca, 'FontSize', 22, 'FontName', 'Times New Roman', 'LineWidth', 1.5);



function [is_c, is_f] = get_time_logic_expanded(t_curr, k, type)

    is_c = false; 
    is_f = false;
    
    if strcmp(type, '30/30')
        
        switch k
            
            case 1
                if (t_curr>12 && t_curr<=42)||(t_curr>66 && t_curr<=96), is_c = true;
                elseif (t_curr>42 && t_curr<=66)||(t_curr>96 && t_curr<=120), is_f = true; end
            
            case 2
                % NEW 1.1% case → using same timing as case 2 (adjust if needed)
                if (t_curr>11 && t_curr<=41)||(t_curr>63 && t_curr<=93), is_c = true;
               elseif (t_curr>41 && t_curr<=63)||(t_curr>93 && t_curr<=115), is_f = true; end

            case 3
                if (t_curr>10 && t_curr<=40)||(t_curr>60 && t_curr<=90), is_c = true;
                elseif (t_curr>40 && t_curr<=60)||(t_curr>90 && t_curr<=110), is_f = true; end
            
            case 4
                if (t_curr>105 && t_curr<=135)||(t_curr>153 && t_curr<=183), is_c = true;
                elseif (t_curr>135 && t_curr<=153)||(t_curr>183 && t_curr<=201), is_f = true; end
            
            case 5
                if (t_curr>560 && t_curr<=590)||(t_curr>606 && t_curr<=636), is_c = true;
                elseif (t_curr>590 && t_curr<=606)||(t_curr>636 && t_curr<=652), is_f = true; end
            
            
        end
        
    elseif strcmp(type, '60/0')

    switch k
        
        case 6   % 1.2%
            if (t_curr > 12 && t_curr <= 72)
                is_c = true;
            elseif (t_curr > 72 && t_curr <= 120)
                is_f = true;
            end

        case 7   % 1.1%
            if (t_curr > 11 && t_curr <= 71)
                is_c = true;
            elseif (t_curr > 71 && t_curr <= 115)
                is_f = true;
            end

        case 8   % 1.0%
            if (t_curr > 10 && t_curr <= 70)
                is_c = true;
            elseif (t_curr > 70 && t_curr <= 110)
                is_f = true;
            end

        case 9   % 0.9%
            if (t_curr > 9 && t_curr <= 69)
                is_c = true;
            elseif (t_curr > 69 && t_curr <= 105)
                is_f = true;
            end

        case 10  % 0.8%
            if (t_curr > 8 && t_curr <= 68)
                is_c = true;
            elseif (t_curr > 68 && t_curr <= 100)
                is_f = true;
            end

    end
        


elseif strcmp(type, '0/60')

    switch k



     case 11   % 1.2%
            if (t_curr > 36 && t_curr <= 96)
                is_c = true;
            elseif (t_curr > 12 && t_curr <= 36) || (t_curr > 96 && t_curr <= 120)
                is_f = true;
            end

        case 12   % 1.1%
            if (t_curr > 127 && t_curr <=187 )
                is_c = true;
            elseif (t_curr > 105 && t_curr <= 127) || (t_curr > 187 && t_curr <= 209)
                is_f = true;
            end

        case 13   % 1.0%
            if (t_curr > 30 && t_curr <= 90)
                is_c = true;
            elseif (t_curr > 10 && t_curr <= 30) || (t_curr > 90 && t_curr <= 110)
                is_f = true;
            end

        case 14   % 0.9%
            if (t_curr > 27 && t_curr <= 87)
                is_c = true;
            elseif (t_curr > 9 && t_curr <= 27) || (t_curr > 87 && t_curr <= 105)
                is_f = true;
            end

        case 15  % 0.8%
            if (t_curr > 100 && t_curr <= 176)
                is_c = true;
            elseif (t_curr > 100 && t_curr <= 116) || (t_curr > 176 && t_curr <= 191)
                is_f = true;
            end
        
    end
    end
end