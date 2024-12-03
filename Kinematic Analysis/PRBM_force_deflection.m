clear; clc; % Define symbolic variables
syms S1 K1 K2 K3 real

% Define parameters
L = 50; % Length of the link in mm
L_flexture = 50;
E = 200;
t = 0.2;
w = 10 ;
theta2_0 = deg2rad(18); % Initial angle in radians (18 degrees)
theta3_0 = 2*pi - theta2_0; % Initial angle for theta3
S_initial = 2 * L * cos(theta2_0); % Initial value of S1 based on theta2_0

% Define the range of S1 within valid bounds
S1_range = linspace(S_initial, S_initial-20, 10);

% Calculate theta2 and theta3 as functions of S1
theta2_formula = acos(S1 / (2 * L));
theta3_formula = 2*pi - theta2_formula;

% Calculate dtheta2/dS1 and dtheta3/dS1 using the given relation
dtheta2_dS1 = 1 / (2 * L * sin(theta2_formula));
dtheta3_dS1 = -1 / (2 * L * sin(theta2_formula));

% Define phi values relative to initial conditions
phi1 = theta2_formula - theta2_0;
phi2 = (theta3_formula - theta3_0) - (theta2_formula - theta2_0);
phi3 = theta3_formula - theta3_0;

% Increase spring constants to match the scale in the target plot

K = (E)*(t*w^3/12)/(L_flexture);
K1 = K;
K2 = K;
K3 = K;

% Substitute constants into the force calculation
Force_expr = phi1 * K1 * dtheta2_dS1 + phi2 * K2 * (dtheta3_dS1 - dtheta2_dS1) + phi3 * K3 * dtheta3_dS1;

% Calculate the initial force when S1 = S_initial
Initial_force = double(subs(Force_expr, S1, S_initial));

% Pre-allocate arrays
theta2_values = zeros(size(S1_range));
theta3_values = zeros(size(S1_range));
force_values = zeros(size(S1_range));

% Loop through each S1 value
for i = 1:length(S1_range)
    S1_value = S1_range(i); % Current S1 value

    % Substitute S1 into the expressions
    Force_sub_S1 = subs(Force_expr, S1, S1_value);
    
    % Evaluate theta2 and theta3
    theta2_values(i) = double(subs(theta2_formula, S1, S1_value));
    theta3_values(i) = double(subs(theta3_formula, S1, S1_value));
    
    % Try evaluating the force and subtract the initial force
    force_values(i) = double(vpa(Force_sub_S1, 10)) - Initial_force;
end

% Plot Force-Deflection graph
figure;
plot(-(S1_range - S_initial), force_values, '--sm'); % Invert the displacement values
xlabel('Deflection S1 (mm)');
ylabel('Force (N)');
title('Force-Deflection Graph');
ylim([0 inf]); % Set y-axis minimum to 0
grid on;

% Display results
disp('Theta2 values (deg):');
disp(rad2deg(theta2_values));
disp('Theta3 values (deg):');
disp(rad2deg(theta3_values));
disp('Force values (N):');
disp(force_values);

% Define data for export
Deflection_S1 = -(S1_range - S_initial); % Inverted displacement values
Force = force_values;               % Scaled forces
Theta2_deg = rad2deg(theta2_values);    % Convert theta2 to degrees
Theta3_deg = rad2deg(theta3_values);    % Convert theta3 to degrees

% Combine data into a table
plot_data = table(Deflection_S1', Force', Theta2_deg', Theta3_deg', ...
    'VariableNames', {'Deflection_S1_mm', 'Force_N', 'Theta2_deg', 'Theta3_deg'});

% Export data to an Excel file
filename = 'PRBM - Force_Deflection_Data.xlsx'; % Specify the file name
writetable(plot_data, filename);

% Notify user
disp(['Plot data has been successfully exported to ', filename]);

% Pre-calculate constants for stress calculations
C = t / 2; % Distance from neutral axis to the surface
I_value = (w * t^3) / 12; % Moment of inertia for rectangular cross-section

% Calculate Moments
M_1 = K1 * phi1;
M_2 = K2 * phi2;
M_3 = K3 * phi3;

% Calculate stresses at each joint
sigma_1 = M_1 * C / I_value;
sigma_2 = M_2 * C / I_value;
sigma_3 = M_3 * C / I_value;

% Evaluate stresses over the range of S1
sigma_1_values = zeros(size(S1_range));
sigma_2_values = zeros(size(S1_range));
sigma_3_values = zeros(size(S1_range));

for i = 1:length(S1_range)
    S1_value = S1_range(i)
    
    % Substitute S1 into the updated phi relations
    phi1_current = subs(phi1, S1, S1_value);
    phi2_current = subs(phi2, S1, S1_value);
    phi3_current = subs(phi3, S1, S1_value);
    
    % Update moments
    M_1_current = double(subs(M_1, S1, S1_value));
    M_2_current = double(subs(M_2, S1, S1_value));
    M_3_current = double(subs(M_3, S1, S1_value));
    
    % Calculate stresses
    sigma_1_values(i) = double(M_1_current * C / I_value);
    sigma_2_values(i) = double(M_2_current * C / I_value);
    sigma_3_values(i) = double(M_3_current * C / I_value);
end

% Plot stress vs. Theta 2
figure;
plot(rad2deg(theta2_values), abs(sigma_1_values), 'r', 'DisplayName', '\sigma_1');
hold on;
plot(rad2deg(theta2_values), sigma_2_values, 'g', 'DisplayName', '\sigma_2');
plot(rad2deg(theta2_values), sigma_3_values, 'b', 'DisplayName', '\sigma_3');
xlabel('Theta2 (degrees)');
ylabel('Stress (Pa)');
title('Stress vs. Theta2');
legend('show');
grid on;

% Notify user
disp('Stress vs. Theta2 plot successfully generated.');





