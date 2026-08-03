% ... [Keep your existing Setup and Global Constants sections] ...

% Initialize arrays to store individual damage components
Dc_values = zeros(1,5);
Df_values = zeros(1,5);
pred_lives_linear = zeros(1,5);
pred_lives_q      = zeros(1,5);
exp_lives         = zeros(1,5);

%% 2. MAIN COMPUTATION LOOP
for k = 1:length(cases)
    eps_a = cases(k).epsilon_a;
    current_file = cases(k).filename;
    B1 = 389501 - (630849 * eps_a) + (257797 * eps_a^2);
    phi = B1 * exp(-Q / (R * T));
    
    if isfile(current_file)
        data = readmatrix(current_file);
        time_col = data(:,1);
        entropy_rate_col = data(:,25);
        Delta_Sf = 0; Dc = 0;
        
        for i = 2:length(time_col)
            t_current = time_col(i);
            dt = time_col(i) - time_col(i-1);
            s_dot = entropy_rate_col(i);
            
            is_creep = false; is_fatigue = false;
            
            % ----- CASE TIME WINDOWS (Logic remains the same) -----
            if k==1
                if (t_current>552 && t_current<=612), is_creep=true;
                elseif (t_current>612 && t_current<=660), is_fatigue=true; end
            elseif k==2
                if (t_current>100 && t_current<=160), is_creep=true;
                elseif (t_current>160 && t_current<=200), is_fatigue=true; end
            elseif k==3
                if (t_current>297 && t_current<=357), is_creep=true;
                elseif (t_current>357 && t_current<=393), is_fatigue=true; end
            elseif k==4
                if (t_current>310 && t_current<=370), is_creep=true;
                elseif (t_current>370 && t_current<=410), is_fatigue=true; end
            elseif k==5
                if (t_current>183 && t_current<=243), is_creep=true;
                elseif (t_current>243 && t_current<=271), is_fatigue=true; end
            end
            
            if is_creep && s_dot > 0
                Dc = Dc + (1/phi)*(s_dot)^(1-n1)*dt;
            elseif is_fatigue
                Delta_Sf = Delta_Sf + s_dot*dt;
            end
        end
        
        % Calculate Fatigue Damage Df using your logarithmic model
        if Delta_Sf >= Sg
            Df = 1;
        else
            Df = (Dfc/log(1-Sc/Sg))*log(1-Delta_Sf/Sg);
        end
        
        % STORE THE INDIVIDUAL VALUES
        Dc_values(k) = Dc;
        Df_values(k) = Df;
        
        % Life Calculation
        pred_lives_linear(k) = 1 / (Dc + Df);
        pred_lives_q(k)      = 1 / (Dc^q + Df^q);
        
    else
        Dc_values(k) = NaN;
        Df_values(k) = NaN;
    end
    exp_lives(k)  = cases(k).exp_life;
end

%% 3. DISPLAY DAMAGE SUMMARY TABLE
T_summary = table((1:5)', [cases.epsilon_a]', Dc_values', Df_values', exp_lives', pred_lives_q', ...
    'VariableNames', {'Case', 'Strain_Amp', 'Creep_Damage_Dc', 'Fatigue_Damage_Df', 'Exp_Life', 'Pred_Life_q'});
disp('--- Creep and Fatigue Damage per Cycle ---');
disp(T_summary);

% ... [Rest of your plotting code] ...
%% ========================================================================
%% 6. PLOT 4: DAMAGE PARTITIONING (D_c vs D_f)
%% ========================================================================
% This plot shows the interaction between creep damage (Y) and fatigue damage (X)
figure('Color','w','Position', [200 200 800 800]); hold on;

% Define colors and markers consistent with your previous logic
colors  = {'r','b','g','m','k'}; 
markers = {'s','d','^','v','p'};
labels = {'\Delta\epsilon = 2.4%', '\Delta\epsilon = 0.8%', '\Delta\epsilon = 1.8%', ...
          '\Delta\epsilon = 2.0%', '\Delta\epsilon = 0.7%'};

% Plot each case
for k = 1:5
    if ~isnan(Dc_values(k))
        loglog(Df_values(k), Dc_values(k), markers{k}, 'MarkerSize', 15, 'LineWidth', 2, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k});
    end
end

% Formatting axes
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 22, 'LineWidth', 1.5, ...
         'FontName', 'Times New Roman', 'TickDir', 'in');

% Determine plot limits
all_d = [Dc_values, Df_values];
d_min = min(all_d(all_d > 0)) * 0.1;
d_max = max(all_d) * 10;
xlim([d_min d_max]); ylim([d_min d_max]);

% Add the 1:1 line (Equal Damage Line)
plot([d_min d_max], [d_min d_max], 'k--', 'LineWidth', 1.5);
text(d_max*0.1, d_max*0.2, 'Creep Dominant \rightarrow', 'FontSize', 18, 'Rotation', 45);
text(d_max*0.2, d_max*0.1, '\leftarrow Fatigue Dominant', 'FontSize', 18, 'Rotation', 45);

% Labels and Title
xlabel('Fatigue Damage per Cycle, D_f', 'FontSize', 26);
ylabel('Creep Damage per Cycle, D_c', 'FontSize', 26);
title('Damage Partitioning (CFI)', 'FontWeight', 'normal');
grid on; box on; axis square;

% Add legend
legend(labels, 'Location', 'southeast', 'FontSize', 16);

hold off;
%% ========================================================================
%% 7. PLOT 5: CREEP-FATIGUE DAMAGE INTERACTION LOCUS
%% ========================================================================
figure('Color','w','Position', [300 300 800 700]); hold on;

% 1. Create the Interaction Locus (Equivalent Damage Line)
% For Linear Damage Summation: Dc + Df = 1
% For NDS (q=0.4): Dc^q + Df^q = 1 -> Dc = (1 - Df^q)^(1/q)
df_line = linspace(0, 1, 1000);
dc_line_linear = 1 - df_line; % Linear Locus
dc_line_nds = (max(0, 1 - df_line.^q)).^(1/q); % NDS Locus

% Plot the Locus curves
plot(df_line, dc_line_linear, 'k--', 'LineWidth', 1.5, 'DisplayName', 'LDS (D_c+D_f=1)');
plot(df_line, dc_line_nds, 'r-', 'LineWidth', 2.5, 'DisplayName', 'NDS (q=0.4) Locus');

% 2. Plot your simulated damage points
% Note: To match the reference figure, these should be cumulative damage 
% (D * N_f) which should theoretically lie on or near the lines.
colors  = {'r','b','g','m','k'}; 
markers = {'o','s','^','d','v'};
labels = {'2.4%','0.8%','1.8%','2.0%','0.7%'};

for k = 1:5
    % Cumulative damage at failure
    Dc_total = Dc_values(k) * pred_lives_q(k);
    Df_total = Df_values(k) * pred_lives_q(k);
    
    plot(Df_total, Dc_total, markers{k}, 'MarkerSize', 12, 'LineWidth', 2, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors{k}, 'DisplayName', labels{k});
end

% 3. Formatting to match Fig. 19 style
set(gca, 'XScale', 'log', 'FontSize', 22, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlabel('Fatigue Damage, D_f', 'FontSize', 26);
ylabel('Creep Damage, D_c', 'FontSize', 26);
xlim([1e-5 2]); ylim([0 1.2]);
grid on; box on;
legend('Location', 'northeastoutside');
title('Creep-Fatigue Interaction Locus', 'FontWeight', 'normal');

hold off;