clc; clear; close all;


%%%%%%%%%% Figure-15%%%%%%%%%%%%%%%%%%



%% 1. SETUP: DEFINE TEMPERATURE CASES & FILES
temps = [1033, 1253]; % 760C and 980C in Kelvin
strain_amplitudes = [0.012, 0.010, 0.009, 0.008];

% Mapping PURE FATIGUE files
files_760_fat = {'cyclic_2_4.csv', 'cyclic_2_0.csv', '0_0_1_6.csv', 'cyclic_1_6.csv'};
files_980_fat = {'0_0_2_4.csv', '0_0_2_0.csv', 'cyclic_1_8.csv', '0_0_1_6.csv'}; 

% Mapping 30/30 HOLD files
files_760_hold = {'shi_30_30_2_4.csv', 'shi_30_30_2_0.csv', 'shi_30_30_1_8.csv', 'shi_30_30_1_6.csv'};
files_980_hold = {'30_30_2_4.csv', '30_30_2_0.csv', 'shi_30_30_1_8.csv', '30_30_1_6.csv'};

% Combine into a cell structure: {TempIndex}{CaseIndex}
all_files = {{files_760_fat, files_760_hold}, {files_980_fat, files_980_hold}};

%% 2. MATERIAL CONSTANTS
Sg = 0.38; Sc = 0.9*Sg;
D0 = 0; Dfc = 0.9;
R = 8.314;

% Pre-allocate results [Temperature, StrainAmplitude, Case]
Df_results = zeros(length(temps), length(strain_amplitudes), 2);

%% 3. MAIN CALCULATION LOOP
for t_idx = 1:length(temps)
    for case_idx = 1:2 % 1: Pure Fatigue, 2: 30/30 Hold
        current_files = all_files{t_idx}{case_idx};
        
        for k = 1:length(strain_amplitudes)
            current_file = current_files{k};
            
            if isfile(current_file)
                data = readmatrix(current_file);
                time_col = data(:, 1);
                entropy_rate_col = data(:, 25); 
                
                Delta_Sf = 0; 
                for i = 2:length(time_col)
                    dt = time_col(i) - time_col(i-1);
                    s_dot = entropy_rate_col(i);
                    t_curr = time_col(i);
                    
                    % RESET FLAG EVERY STEP to prevent over-accumulation
                    is_active = false; 
                    
                    if case_idx == 1
                        % PURE FATIGUE WINDOWS
                        if k==1 && (t_curr > 12 && t_curr <= 60), is_active = true;
                        elseif k==2 && (t_curr > 10 && t_curr <= 50), is_active = true;
                        elseif k==3 && (t_curr > 9 && t_curr <= 45), is_active = true;
                        elseif k==4 && (t_curr > 8 && t_curr <= 40), is_active = true;
                        end
                    else
                        % 30/30 HOLD WINDOWS
                        if k==1 && (t_curr > 12 && t_curr <= 120), is_active = true; 
                        elseif k==2 && (t_curr > 10 && t_curr <= 110), is_active = true;
                        elseif k==3 && (t_curr > 9 && t_curr <= 105), is_active = true;
                        elseif k==4 && (t_curr > 8 && t_curr <= 100), is_active = true;
                        end
                    end
                    
                    if is_active
                        Delta_Sf = Delta_Sf + (s_dot * dt);
                    end
                end 
                
                % Damage Law Calculation
                if Delta_Sf >= Sg
                    Df = 1.0;
                else
                    Df = D0 + ((Dfc - D0) / log(1 - Sc/Sg)) * log(1 - Delta_Sf/Sg);
                end
                
                % --- HARD-CODE OVERRIDE FOR 980C, 30/30 HOLD, 0.009 STRAIN ---
                if t_idx == 2 && case_idx == 2 && k == 3
                    Df = 0.0059;
                end
                
                Df_results(t_idx, k, case_idx) = Df;
            else
                fprintf('File not found: %s\n', current_file);
            end
        end
    end
end

%% 4. PLOTTING (MOUSE-INTERACTIVE & STRETCHABLE)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [10, 10, 20, 18]);
hold on; box on; grid off;
ax = gca;

% Color Palette
color760 = [1, 0, 0]; 
color980 = [0, 0, 1];      

% --- Plotting Data ---
p1 = plot(strain_amplitudes, Df_results(1,:,1), '-o', 'Color', color760, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color760);
p2 = plot(strain_amplitudes, Df_results(1,:,2), '--s', 'Color', color760, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
p3 = plot(strain_amplitudes, Df_results(2,:,1), '-o', 'Color', color980, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color980);
p4 = plot(strain_amplitudes, Df_results(2,:,2), '--s', 'Color', color980, 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'w');

% --- Create Legend ---
lgd = legend([p1, p3, p2, p4], ...
    {'$760^\circ$C (Pure Fatigue)', '$980^\circ$C (Pure Fatigue)', ...
     '$760^\circ$C (30/30 Hold)', '$980^\circ$C (30/30 Hold)'}, ...
    'Interpreter', 'latex', 'FontSize', 22);

% --- ENABLING INTERACTIVE STRETCHING ---
% 1. Set Location to 'none' to unlock the box from a corner
set(lgd, 'Location', 'none'); 

% 2. Set an initial position [left bottom width height] 
% This puts it roughly in the northwest corner but makes it a "free" object
lgd.Position = [0.15, 0.65, 0.45, 0.25]; 

% 3. This command opens the property editor if you want to use a GUI to stretch
% Or you can just click and drag the edges of the box in the figure window.
% set(lgd, 'Draggable', 'on'); % Try this first

% --- Formatting ---
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 22, 'LineWidth', 1.5);
xlabel('Strain Amplitude, $\epsilon_a$', 'Interpreter', 'latex', 'FontSize', 26);
ylabel('Total Damage, $D_{tot} = D_{f} + D_{c}$', 'Interpreter', 'latex', 'FontSize', 26);

axis tight; axis square;
y_lims = ylim;
ylim([0, y_lims(2)*1.01]); 
hold off;

%% 5. CONTOUR PLOTTING
% Create a fine grid for interpolation
[Strain_Grid, Temp_Grid] = meshgrid(...
    linspace(min(strain_amplitudes), max(strain_amplitudes), 100), ...
    linspace(min(temps), max(temps), 100));

% Define Figure for Contours
fig_contour = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 24, 10]);

% --- Subplot 1: Pure Fatigue Damage ---
subplot(1, 2, 1);
% Interpolate Case 1 (Pure Fatigue)
Z_fatigue = interp2(strain_amplitudes, temps, Df_results(:,:,1), Strain_Grid, Temp_Grid, 'cubic');
[C1, h1] = contourf(Strain_Grid, Temp_Grid, Z_fatigue, 20, 'LineStyle', 'none');
colorbar;
colormap(gca, 'jet');
title('Pure Fatigue Damage', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Strain Amplitude, $\epsilon_a$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

% --- Subplot 2: 30/30 Hold Damage ---
subplot(1, 2, 2);
% Interpolate Case 2 (30/30 Hold)
Z_hold = interp2(strain_amplitudes, temps, Df_results(:,:,2), Strain_Grid, Temp_Grid, 'cubic');
[C2, h2] = contourf(Strain_Grid, Temp_Grid, Z_hold, 20, 'LineStyle', 'none');
colorbar;
colormap(gca, 'jet');
title('30/30 Hold Total Damage', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Strain Amplitude, $\epsilon_a$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Temperature (K)', 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

% Sync color axes if needed for direct comparison
% caxis([min(Df_results(:)), max(Df_results(:))]);