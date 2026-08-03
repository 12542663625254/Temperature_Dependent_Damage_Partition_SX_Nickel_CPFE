clc; clear; close all;

%% 1. SETUP: DEFINE TEMPERATURE CASES & FILES
temps = [1033, 1253]; % 760C and 980C in Kelvin
strain_amplitudes = [0.012, 0.010, 0.009, 0.008];

% Mapping CSV files for each temperature
files_760 = {'cyclic_2_4.csv', 'cyclic_2_0.csv', '0_0_1_6.csv', 'cyclic_1_6.csv'};
files_980 = {'0_0_2_4.csv', '0_0_2_0.csv', 'cyclic_1_8.csv', '0_0_1_6.csv'}; 
all_files = {files_760, files_980};

%% 2. MATERIAL CONSTANTS
Sg = 0.38; Sc = 0.9*Sg;
D0 = 0; Dfc = 0.9;
n1 = 0.6; Q = 6.97e-19; R = 8.314;

% Pre-allocate results for plotting
Df_results = zeros(length(temps), length(strain_amplitudes));

%% 3. MAIN CALCULATION LOOP
for t_idx = 1:length(temps)
    T = temps(t_idx);
    current_temp_files = all_files{t_idx};
    
    for k = 1:length(strain_amplitudes)
        eps_a = strain_amplitudes(k);
        current_file = current_temp_files{k};
        
        if isfile(current_file)
            data = readmatrix(current_file);
            time_col = data(:, 1);
            entropy_rate_col = data(:, 25); 
            
            Delta_Sf = 0; 
            
            for i = 2:length(time_col)
                dt = time_col(i) - time_col(i-1);
                s_dot = entropy_rate_col(i);
                t_curr = time_col(i);
                
                % Fatigue window logic
                is_fatigue = false;
                if k==1 && (t_curr > 12 && t_curr <= 60), is_fatigue = true;
                elseif k==2 && (t_curr > 10 && t_curr <= 50), is_fatigue = true;
                elseif k==3 && (t_curr > 9 && t_curr <= 45), is_fatigue = true;
                elseif k==4 && (t_curr > 8 && t_curr <= 40), is_fatigue = true;
                end
                
                if is_fatigue
                    Delta_Sf = Delta_Sf + (s_dot * dt);
                end
            end
            
            % Logarithmic Fatigue Damage law
            if Delta_Sf >= Sg
                Df = 1.0;
            else
                Df = D0 + ((Dfc - D0) / log(1 - Sc/Sg)) * log(1 - Delta_Sf/Sg);
            end
            
            Df_results(t_idx, k) = Df;
        else
            fprintf('Warning: File %s not found.\n', current_file);
        end
    end
end

%% 4. PLOTTING FATIGUE DAMAGE (PROFESSIONAL QUALITY - ENHANCED)
fig = figure('Color', 'w');
fig.Units = 'centimeters';
fig.Position = [10, 10, 16, 14]; % Slightly larger for higher visibility

hold on;
box on; 
grid off;
ax = gca;

% Professional color palette
color1 = [0.8500, 0.3250, 0.0980]; % Deep Red
color2 = [0, 0.4470, 0.7410];      % Deep Blue

% --- Plotting with Heavy LineWeights and Markers ---
% Plot 760C - Solid line, bold markers
p1 = plot(strain_amplitudes, Df_results(1,:), 'r', ...
    'Color', color1, ...
    'LineWidth', 3, ...        % Increased line width
    'MarkerSize', 11, ...      % Increased marker size
    'MarkerFaceColor', color1, ...
    'MarkerEdgeColor', 'k');

% Plot 980C - Dashed line, bold open markers
p2 = plot(strain_amplitudes, Df_results(2,:), 'b', ...
    'Color', color2, ...
    'LineWidth', 3, ...        % Increased line width
    'MarkerSize', 12, ...      % Increased marker size
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', color2, ...
    'MarkerSize', 12);

% --- Axis Formatting ---
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 20);
ax.LineWidth = 1.5; % Bold axis frame
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [0.02, 0.02]; % Longer ticks for a technical look

% Labeling
xlabel('Strain Amplitude, $\epsilon_a$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('Fatigue Damage, $D_f$', 'Interpreter', 'latex', 'FontSize', 24);

% --- Legend Styling (Boxed) ---
lgd = legend([p1, p2], {'$T = 760$ $^\circ$C', '$T = 980$ $^\circ$C'}, ...
    'Location', 'best', ...
    'Interpreter', 'latex', ...
    'FontSize', 20, ...
    'EdgeColor', [0.2 0.2 0.2], ... % Dark grey border
    'LineWidth', 1);               % Box border thickness

% --- Final Polish ---
axis tight;
axis square; % Makes the plot perfectly square
y_limits = ylim;
% Adding 15% padding to the top for better breathing room
ylim([y_limits(1), y_limits(2)*1.15]); 

hold off;

% Export command for publication-ready PDF
% exportgraphics(fig, 'Fatigue_Damage_Final.pdf', 'ContentType', 'vector', 'Colorspace', 'rgb');