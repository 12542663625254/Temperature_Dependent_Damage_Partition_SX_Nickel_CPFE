% 1. Load the data
data = readmatrix('Fatigue.csv');

% 2. Extract specific columns
% Shift x-values by subtracting 103 so the plot "starts" at 0
x = data(:, 1) - 103.5; 
y = data(:, 36); 

% 3. Create the plot
figure('Color', 'w');
plot(x, y, 'LineWidth', 2, 'Color', [0, 0.447, 0.741]); 

% Adjust x-limits to match the shift (265 - 103 = 162)
xlim([0 162]);
hold on; 

% 4. Add horizontal midline
y_limits = ylim;
y_mid = (y_limits(1) + y_limits(2)) / 2;
yline(y_mid, '--r', 'LineWidth', 2);

% 5. Add professional formatting with Times New Roman
% Set font properties for the axes (tick marks)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22);

% Set labels and title with Times New Roman
xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 26);
ylabel('Stress (MPa)', 'FontName', 'Times New Roman', 'FontSize', 26);


grid off;
box on;
hold off;