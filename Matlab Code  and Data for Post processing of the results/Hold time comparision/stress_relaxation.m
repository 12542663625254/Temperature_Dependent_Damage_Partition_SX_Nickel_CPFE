clc; clear; close all;

%% ============================================================
%% 1. FILE NAMES
%% ============================================================
file_60  = 'Fatigue_60_0_2_0.csv';
file_240 = '1_0_strain.csv';   % <-- Replace if needed

%% ============================================================
%% 2. TIME WINDOWS
%% ============================================================
t_start     = 10;   % Hold start time
t_end_60    = 70;   % 60 s hold end
t_end_240   = 40;   % Adjust based on your second dataset

%% ============================================================
%% 3. CREATE FIGURE
%% ============================================================
figure('Color', 'w', 'Position', [100 100 950 700]);
hold on;

%% ============================================================
%% 4. PROCESS 60s DATASET
%% ============================================================
if isfile(file_60)
    
    data60 = readmatrix(file_60);
    time60_raw   = data60(:,1);
    stress60_raw = data60(:,36);

    idx60 = (time60_raw >= t_start) & (time60_raw <= t_end_60);
    time60   = time60_raw(idx60) - t_start;
    stress60 = stress60_raw(idx60);

    % Plot curve
    plot(time60, stress60, 'r-', 'LineWidth', 2.5);

    % Stress drop
    delta_sigma60 = stress60(1) - stress60(end);

    % Start & End markers
    plot(time60(1),  stress60(1),  'ro', ...
        'MarkerSize',10,'MarkerFaceColor','r');
    plot(time60(end), stress60(end), 'rs', ...
        'MarkerSize',10,'MarkerFaceColor','r');

    % Vertical drop line
    plot([time60(end) time60(end)], ...
         [stress60(1) stress60(end)], ...
         'r--', 'LineWidth',1.8);

    % Annotation
    text(time60(end)+5, mean([stress60(1) stress60(end)]), ...
        sprintf('$\\Delta\\sigma_{60} = %.1f$ MPa', delta_sigma60), ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'Color','r', ...
        'FontWeight','bold');

else
    warning('60s file not found');
end

%% ============================================================
%% 5. PROCESS SECOND DATASET (240s or Other)
%% ============================================================
if isfile(file_240)
    
    data240 = readmatrix(file_240);
    time240_raw   = data240(:,1);
    stress240_raw = data240(:,36);

    idx240 = (time240_raw >= t_start) & (time240_raw <= t_end_240);
    time240   = time240_raw(idx240) - t_start;
    stress240 = stress240_raw(idx240);

    % Plot curve
    plot(time240, stress240, 'b-', 'LineWidth', 2.5);

    % Stress drop
    delta_sigma240 = stress240(1) - stress240(end);

    % Start & End markers
    plot(time240(1),  stress240(1),  'bo', ...
        'MarkerSize',10,'MarkerFaceColor','b');
    plot(time240(end), stress240(end), 'bs', ...
        'MarkerSize',10,'MarkerFaceColor','b');

    % Vertical drop line
    plot([time240(end) time240(end)], ...
         [stress240(1) stress240(end)], ...
         'b--', 'LineWidth',1.8);

    % Annotation
    text(time240(end)+5, mean([stress240(1) stress240(end)]), ...
        sprintf('$\\Delta\\sigma_{240} = %.1f$ MPa', delta_sigma240), ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'Color','b', ...
        'FontWeight','bold');

else
    warning('Second dataset file not found');
end

%% ============================================================
%% 6. FORMATTING
%% ============================================================
set(gca, 'FontSize', 22, ...
         'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, ...
         'TickDir', 'out');

grid on; box on;

xlabel('Hold Time, $t_{hold}$ (s)', ...
    'Interpreter', 'latex', ...
    'FontSize', 26, ...
    'FontWeight','bold');

ylabel('Stress, $\sigma$ (MPa)', ...
    'Interpreter', 'latex', ...
    'FontSize', 26, ...
    'FontWeight','bold');

title('Stress Relaxation Comparison Under Different Hold Durations', ...
    'FontSize', 24, ...
    'FontWeight','bold');

xlim([0 250]);
ylim([500 1100]);

legend({'60 s Hold','Second Dataset'}, ...
    'FontSize',20, ...
    'Location','northeast');

hold off;
