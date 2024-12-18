clear; clc; % Clear workspace and command window

% Define symbolic variables
syms X K1 K2 K3 real

% Define parameters
L = 112.5; % Length of the link in mm
L_flexture = 50; % Length of the flexure
E = 200; % Young's modulus (assumed units consistent with inputs)
t = 0.4; % Thickness in mm
w = 30; % Width in mm
theta_20 = deg2rad(30); % Initial angle theta_2 in radians (30 degrees)
theta_30 = 2*pi - theta_20; % Initial angle theta_3
X_initial = 2 * L * cos(theta_20); % Initial value of X based on theta_20

% Define the range of X within valid bounds
X_range = linspace(X_initial, X_initial - 20, 10);

% Calculate theta_2 and theta_3 as functions of X
theta_2_formula = acos(X / (2 * L));
theta_3_formula = 2*pi - theta_2_formula;

% Calculate d(theta_2)/dX and d(theta_3)/dX:
dtheta_2_dX = 1 / (2 * L * sin(theta_2_formula));
dtheta_3_dX = -1 / (2 * L * sin(theta_2_formula));

% Define phi values relative to initial conditions
phi1 = theta_2_formula - theta_20;
phi2 = (theta_3_formula - theta_30) - (theta_2_formula - theta_20);
phi3 = theta_3_formula - theta_30;

% Calculate spring constants
K = (E) * (t * w^3 / 12) / (L_flexture); % Simplified flexure stiffness equation
K1 = K;
K2 = K;
K3 = K;

% Force expression
Force_expr = phi1 * K1 * dtheta_2_dX + phi2 * K2 * (dtheta_3_dX - dtheta_2_dX) + phi3 * K3 * dtheta_3_dX;

% Calculate the initial force when X = X_initial
Initial_force = double(subs(Force_expr, X, X_initial));

% Pre-allocate arrays
theta_2_values = zeros(size(X_range));
theta_3_values = zeros(size(X_range));
force_values = zeros(size(X_range));

% Loop through each X value
for i = 1:length(X_range)
    X_value = X_range(i); % Current X value

    % Substitute X into the expressions
    Force_sub_X = subs(Force_expr, X, X_value);

    % Evaluate theta_2 and theta_3
    theta_2_values(i) = double(subs(theta_2_formula, X, X_value));
    theta_3_values(i) = double(subs(theta_3_formula, X, X_value));

    % Evaluate force and subtract the initial force
    force_values(i) = double(vpa(Force_sub_X, 10)) - Initial_force;
end

% Plot Force-Deflection graph
figure;
plot(-(X_range - X_initial), force_values, LineWidth=3); 
xlabel('Deflection X (mm)');
ylabel('Force (N)');
title('Force-Deflection Graph');
ylim([0 inf]);
grid on;

% Display results
disp('Theta_2 values (deg):');
disp(rad2deg(theta_2_values));
disp('Theta_3 values (deg):');
disp(rad2deg(theta_3_values));
disp('Force values (N):');
disp(force_values);

% Define data for export
Deflection_X = -(X_range - X_initial); % Inverted displacement values
Force = force_values;               % Scaled forces
Theta2_deg = rad2deg(theta_2_values);    % Convert theta_2 to degrees
Theta3_deg = rad2deg(theta_3_values);    % Convert theta_3 to degrees

% Combine data into a table
plot_data = table(Deflection_X', Force', Theta2_deg', Theta3_deg', ...
    'VariableNames', {'Deflection_X_mm', 'Force_N', 'Theta2_deg', 'Theta3_deg'});

% Export data to an Excel file
filename = 'Data.csv'; % Specify the file name
writetable(plot_data, filename);
