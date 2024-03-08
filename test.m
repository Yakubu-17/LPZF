% Create the main figure
figure;

% Plot something in the main figure (e.g., a sinusoidal function)
x = linspace(0, 2*pi, 100);
y = sin(x);
plot(x, y, 'b-', 'LineWidth', 2);
title('Main Figure');

% Manually create a smaller axes within the main figure and adjust its position
smaller_axes = axes('Position', [0.6, 0.2, 0.3, 0.3]); % [left, bottom, width, height]
x_sub = linspace(-2, 2, 50);
y_sub = x_sub.^2;
plot(x_sub, y_sub, 'r--', 'LineWidth', 1.5);
title('Smaller Axes');

% Add labels and other decorations to the smaller axes
xlabel('X-axis (Smaller Axes)');
ylabel('Y-axis (Smaller Axes)');
grid on;

% Adjust the layout to avoid overlapping
sgtitle('Main Figure with Smaller Axes');
