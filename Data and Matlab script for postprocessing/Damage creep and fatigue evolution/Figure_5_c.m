clc; clear; close all;

%%%%%% Figure-5(c)%%%%%%%%

%% 1. LOAD DATA  (Case 1: 1.2% strain amplitude, 30/30 hold)
filename  = '30_30_2_4.csv';
epsilon_a = 0.012;

% --- TEP model constants (from your damage script) ---
n1 = 0.6; D0 = 0;
Sg_fixed = 0.38; Sc_ratio_fixed = 0.9; Dfc_fixed = 0.9;
Q = 6.97e-19; T = 1033; R_univ = 8.314;

if isfile(filename)
    data      = readmatrix(filename);
    t         = data(:, 1);    % Time
    s_dot_vec = data(:, 25);   % Entropy generation rate (Column 25)
else
    % Fallback synthetic data so the code still runs for previewing layout
    t         = (0:0.1:120)';
    s_dot_vec = 1e-4 * (1 + 0.5*sin(2*pi*t/60)) .* (t > 0);
end

% --- Force time to start at 0 so data aligns with the partition points ---
% (no-op if the file already starts at 0; corrects a 103.2 s offset if present)
t = t - t(1);

%% 2. PARTITION BOUNDARIES (6 cycle phases, as in the reference figure)
% (1) loading 0-17 | (2) tensile hold 17-47 | (3) unloading 47-54.2 |
% (4) compressive loading 54.2-71 | (5) compressive hold 71-101 | (6) reloading 101-108.2
t_points_absolute = [0, 17, 47, 54.2, 71, 101, 108.2];   % 7 boundaries -> 6 partitions
n_portions        = numel(t_points_absolute) - 1;
t_end             = t_points_absolute(end);              % end of partition 6 = 108.2 s
portion_names     = {'Loading','Tensile hold','Unloading', ...
                     'Compressive loading','Compressive hold','Reloading'};

%% 3. COMPUTE CUMULATIVE Dc & Df + PER-PORTION BREAKDOWN (TEP model)
B1  = 2.0428e3 + (-2.1056e4 * epsilon_a) + (-1.2398e7 * epsilon_a^2);
phi = B1 * exp(-Q / (R_univ * T));

N = length(t);
Dc_series = zeros(N,1);   % cumulative creep damage up to t(i)
Df_series = zeros(N,1);   % fatigue damage from cumulative entropy up to t(i)

Dc_portion = zeros(1, n_portions);   % creep damage accrued in each portion
Df_portion = zeros(1, n_portions);   % fatigue damage accrued in each portion

Delta_Sf = 0; C_sum = 0;
val_in       = max(1e-10, 1 - Delta_Sf / Sg_fixed);
Df_series(1) = D0 + ((Dfc_fixed - D0) / log(1 - Sc_ratio_fixed)) * log(val_in);

for i = 2:N
    dt = t(i) - t(i-1); s_dot = s_dot_vec(i);
    [is_c, is_f] = get_time_logic(t(i));
    p = get_portion(t(i), t_points_absolute);   % which of the 6 portions (0 if outside)

    dDc = 0;
    if is_c && s_dot > 0
        dDc   = (1/phi) * (s_dot^(1-n1) * dt);
        C_sum = C_sum + dDc;
    elseif is_f
        Delta_Sf = Delta_Sf + (s_dot * dt);
    end

    Dc_series(i) = C_sum;
    val_in       = max(1e-10, 1 - Delta_Sf / Sg_fixed);
    Df_series(i) = D0 + ((Dfc_fixed - D0) / log(1 - Sc_ratio_fixed)) * log(val_in);

    dDf = Df_series(i) - Df_series(i-1);   % fatigue increment this step

    % --- tally into the active portion ---
    if p >= 1 && p <= n_portions
        Dc_portion(p) = Dc_portion(p) + dDc;
        Df_portion(p) = Df_portion(p) + dDf;
    end
end

%% 4. REPORT
fprintf('\nFinal creep damage   Dc = %.4e\n',   Dc_series(end));
fprintf('Final fatigue damage Df = %.4e\n\n', Df_series(end));

fprintf('%-4s %-22s %14s %14s\n', 'No.', 'Portion', 'Dc (creep)', 'Df (fatigue)');
fprintf('%s\n', repmat('-',1,56));
for p = 1:n_portions
    fprintf('%-4d %-22s %14.4e %14.4e\n', p, portion_names{p}, Dc_portion(p), Df_portion(p));
end
fprintf('%s\n', repmat('-',1,56));
fprintf('%-4s %-22s %14.4e %14.4e\n\n', '', 'TOTAL', sum(Dc_portion), sum(Df_portion));

%% 5. PLOT SETTINGS
t_start   = 0;
arrow_gap = 1.0;

% --- Mask the CURVE data to end of partition 6 (axis stays full width) ---
data_mask = (t >= t_start) & (t <= t_end);

%% 6. COMBINED FIGURE: Dc (left axis) & Df (right axis) vs TIME
figure('Name', 'Creep & Fatigue Damage', 'Color', 'w', 'Position', [200, 80, 900, 760]);
hold on; box on; axis square;

dc_max = max(Dc_series); if dc_max == 0, dc_max = 1; end
df_max = max(Df_series); if df_max == 0, df_max = 1; end

% --- Left axis: Creep damage (data only up to partition 6) ---
yyaxis left
hC = plot(t(data_mask), Dc_series(data_mask), 'r-', 'LineWidth', 2.5);
ylim([0, dc_max * 1.35]);
ylabel('Creep damage, $D_c$', 'FontSize', 26, 'Interpreter', 'latex');
set(gca, 'YColor', 'r');

% --- Partition markers (drawn in LEFT-axis coordinates) ---
plot_common_elements(t_points_absolute, dc_max*1.10, dc_max*1.22, arrow_gap);

% --- Right axis: Fatigue damage (data only up to partition 6) ---
yyaxis right
hF = plot(t(data_mask), Df_series(data_mask), 'b-', 'LineWidth', 2.5);
ylim([0, df_max * 1.35]);
ylabel('Fatigue damage, $D_f$', 'FontSize', 26, 'Interpreter', 'latex');
set(gca, 'YColor', 'b');

% --- X-axis + shared formatting (FULL window, like the reference figure) ---
xlim([t_start, t_start + 130]);
ticks_abs = t_start : 20 : (t_start + 120);
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman', 'LineWidth', 1.5, ...
    'XTick', ticks_abs, 'XTickLabel', 0:20:120, 'TickDir', 'in');
xlabel('Time (s)', 'FontSize', 26, 'Interpreter', 'latex');

legend([hC, hF], {'$D_c$ (Creep)', '$D_f$ (Fatigue)'}, ...
    'Interpreter', 'latex', 'FontSize', 20, 'Location', 'northwest');

% --- Annotation: hold configuration + strain amplitude ---
ann_str = {'Hold: 30/30 s', ...
           sprintf('$\\epsilon_a = %.1f\\%%$', epsilon_a*100)};
annotation('textbox', [0.58, 0.16, 0.30, 0.12], 'String', ann_str, ...
    'Interpreter', 'latex', 'FontSize', 22, 'FontName', 'Times New Roman', ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

%% --- HELPER FUNCTIONS ---
function p = get_portion(t_curr, t_abs)
    % Returns the portion index (1..6) that t_curr falls into; 0 if outside.
    p = 0;
    for j = 1:(numel(t_abs)-1)
        if t_curr > t_abs(j) && t_curr <= t_abs(j+1)
            p = j; return;
        end
    end
end

function plot_common_elements(t_abs, arrow_y, label_y, gap)
    for i = 1:length(t_abs)
        xline(t_abs(i), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        if i < length(t_abs)
            x1  = t_abs(i) + gap; x2 = t_abs(i+1) - gap;
            mid = (t_abs(i) + t_abs(i+1)) / 2;
            plot([x1, x2], [arrow_y, arrow_y], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(x1, arrow_y, 'k<', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
            plot(x2, arrow_y, 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
            plot(mid, label_y, 'ko', 'MarkerSize', 20, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
            text(mid, label_y, num2str(i), 'FontSize', 12, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end
    end
end

%% HELPER: Time logic aligned to the 6 cycle phases / partition boundaries
% Creep accrues during the two HOLDS  -> portions 2 (17-47 s) & 5 (71-101 s)
% Fatigue accrues during the RAMPS    -> portions 1, 3, 4, 6
function [is_c, is_f] = get_time_logic(t_curr)
    is_c = false; is_f = false;
    if (t_curr > 17  && t_curr <= 47) || (t_curr > 71 && t_curr <= 101)
        is_c = true;                                   % tensile / compressive holds
    elseif (t_curr > 0   && t_curr <= 17)  || ...
           (t_curr > 47  && t_curr <= 71)  || ...
           (t_curr > 101 && t_curr <= 108.2)
        is_f = true;                                   % loading / unloading / reloading ramps
    end
end
