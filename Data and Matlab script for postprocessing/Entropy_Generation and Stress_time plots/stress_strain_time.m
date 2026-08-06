clc; clear; close all;

%%%%%%% Figure-5 (a) and (b) of the present work%%%%%%%%%%%%


%% 1. LOAD DATA
filename = 'stress_strain_30_30_2_4.csv'; 
if isfile(filename)
    data = readmatrix(filename);
    x_data = data(:, 1);    
    y_stress = data(:, 36); 
    y_strain = data(:, 33); % Strain from column 33
    y_entropy = data(:, 25); 
else
    x_data = linspace(103.2, 233.2, 1000);
    y_stress = 1000 * sin((x_data-103.2)/20);
    y_strain = 0.02 * sin((x_data-103.2)/20 - 0.5); 
    y_entropy = 0.05 * (x_data-103.2).^1.2; 
end

%% 2. SET SHARED PARAMETERS
t_start = 103.2;
t_points_relative = [0, 17, 47, 54.2, 71, 101, 108.2]; 
t_points_absolute = t_start + t_points_relative; 
arrow_gap = 1.5;

%% 3. FIGURE 1: STRESS & STRAIN vs TIME
figure('Name', 'Stress and Strain Plot', 'Color', 'w', 'Position', [100, 100, 900, 700]);
hold on;
axis square;

% --- PRE-CALCULATE DATA MODIFICATIONS ---
% Create a copy of strain and set negative values before t=108 to NaN
y_strain_modified = y_strain;
mask_to_hide = (x_data < 108) & (y_strain < 0);
y_strain_modified(mask_to_hide) = NaN;

% --- PRE-CALCULATE LIMITS TO ALIGN ZEROS ---
stress_max = max(y_stress);
stress_min = min(y_stress);
strain_max = max(y_strain, [], 'omitnan');
strain_min = min(y_strain, [], 'omitnan');

% Aligning zeros: Ratio of Top/Bottom must be equal for both axes
s_top = stress_max * 1.6; 
s_bot = stress_min * 1.2;
ratio = s_top / s_bot;
sn_top = strain_max * 1.6;
sn_bot = sn_top / ratio; 

% --- DEFINE OPPOSITE COLOR FOR STRESS ---
stress_color = [1, 0, 1]; % Deep Burnt Orange (Complementary to Blue)

% Left Axis: Stress
yyaxis left
plot(x_data, y_stress, '-', 'Color', stress_color, 'LineWidth', 2.5); 
ylabel('Stress (MPa)', 'FontSize', 26, 'Interpreter', 'latex', 'Color', stress_color);
set(gca, 'YColor', stress_color); % Sets axis, ticks, and labels to the new color
ylim([s_bot, s_top]);

% --- Dynamic positioning for arrows/labels ---
arrow_y_stress = stress_max + 0.1 * (stress_max - stress_min);
label_y_stress = stress_max + 0.2 * (stress_max - stress_min);

% Right Axis: Strain
yyaxis right
plot(x_data, y_strain_modified, 'b-', 'LineWidth', 2.0); 
ylabel('Strain', 'FontSize', 26, 'Interpreter', 'latex', 'Color', 'b');
set(gca, 'YColor', 'b'); 
ylim([sn_bot, sn_top]); 

% ⚠️ CRITICAL FIX: switch back to LEFT axis before numbering
yyaxis left
plot_common_elements(t_points_absolute, arrow_y_stress, label_y_stress, arrow_gap);
yline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
format_plot(t_start, []);
%% 4. FIGURE 2: ENTROPY & STRAIN vs TIME
figure('Name', 'Entropy and Strain Plot', 'Color', 'w', 'Position', [950, 100, 900, 700]);
hold on;
axis square;

% --- PRE-CALCULATE LIMITS TO ALIGN ZEROS ---
ent_max = max(y_entropy);
ent_min = min(y_entropy);
strain_max = max(y_strain, [], 'omitnan');
strain_min = min(y_strain, [], 'omitnan');

% --- ZOOM OUT ADJUSTMENT ---
e_top = ent_max * 1.6; 
e_bot = -1.2 * ent_max; % Deeper floor to show negative strain cycles
ratio_ent = e_top / e_bot;
sn_top_ent = strain_max * 1.6;
sn_bot_ent = sn_top_ent / ratio_ent; % Proportional scaling for right axis

% Left Axis: Entropy
yyaxis left
plot(x_data, y_entropy, 'r-', 'LineWidth', 2.5); 
ylabel('Entropy Gen. Rate ($\dot{S}$) (W/m$^3\cdot$K)', 'Interpreter', 'latex', 'FontSize', 26, 'Color', 'r');

% --- ADDED: DEEP GREEN HORIZONTAL LINE AT 2e-4 ---
yline(2e-4, '--', 'Color', [0, 0.45, 0], 'LineWidth', 2.5, 'HandleVisibility', 'off'); 

set(gca, 'YColor', 'r'); 
ylim([e_bot, e_top]);

% --- Dynamic positioning for arrows (based on Entropy scale) ---
arrow_y_ent = ent_max + 0.1 * (ent_max - ent_min);
label_y_ent = ent_max + 0.3 * (ent_max - ent_min);

% Right Axis: Strain
yyaxis right
y_strain_mod_ent = y_strain;
mask_hide_ent = (x_data < 108) & (y_strain < 0);
y_strain_mod_ent(mask_hide_ent) = NaN;
plot(x_data, y_strain_mod_ent, 'b-', 'LineWidth', 2.0); 
ylabel('Strain', 'FontSize', 26, 'Interpreter', 'latex', 'Color', 'b');
set(gca, 'YColor', 'b'); 
ylim([sn_bot_ent, sn_top_ent]); 

% ⚠️ CRITICAL: switch back to LEFT axis before numbering
yyaxis left
plot_common_elements(t_points_absolute, arrow_y_ent, label_y_ent, arrow_gap);

% Common reference line
yline(0, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off'); 
format_plot(t_start, []);
%% ==========================================================
% LOCAL FUNCTIONS
%% ==========================================================

function plot_common_elements(t_abs, arrow_y, label_y, gap)
    for i = 1:length(t_abs)
        xline(t_abs(i), 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        if i < length(t_abs)
            x1 = t_abs(i) + gap; 
            x2 = t_abs(i+1) - gap;
            mid = (t_abs(i) + t_abs(i+1)) / 2;

            plot([x1, x2], [arrow_y, arrow_y], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(x1, arrow_y, 'k<', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
            plot(x2, arrow_y, 'k>', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');

            plot(mid, label_y, 'ko', 'MarkerSize', 20, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
            text(mid, label_y, num2str(i), ...
                'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
end

function format_plot(t_start, y_lims)
    grid off; box on; axis square;

    xlim([t_start, t_start + 130]);

    if ~isempty(y_lims)
        ylim(y_lims);
    end

    ticks_abs = t_start : 20 : (t_start + 120);

    set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 1.5, ...
             'XTick', ticks_abs, ...
             'XTickLabel', 0:20:120, ...
             'TickDir', 'in');

    xlabel('Time (s)', 'FontSize', 26, 'Interpreter', 'latex');
end