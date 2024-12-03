% Clear workspace and close all figures
clc;
clear;
close all;

%% Symbolic Variables
% Define symbolic variables for angles (theta2, theta3), lengths (L2, L3),
% spring constants (K1, K2, K3), and initial angles (theta2_0, theta3_0)
syms theta2 theta2_0 L2 L3 K1 K2 K3 theta3_0 real

%% Input Force Expression
% The input force F_in is calculated as a function of theta2 and other parameters
theta3 = asin((L2 / L3) * sin(theta2));  % Relationship between theta2 and theta3

% Define the input force equation (symbolic)
F_in = ...
    K1 * (theta2 - theta2_0) * ...
        (-cos(asin((L2 / L3) * sin(theta2))) / (L2 * sin(theta3 - theta2))) + ...
    K2 * (((theta2 - theta2_0) - (theta3 - theta3_0)) * ...
        (-cos(theta3) / (L2 * sin(theta3 - theta2))) + ...
        cos(theta2) / (L3 * sin(theta3 - theta2))) + ...
    K3 * ((theta3 - theta3_0) * -cos(theta2) / (L3 * sin(theta3)));

% Substitute the expression for theta3 into F_in
F_in = subs(F_in, theta3, asin((L2 / L3) * sin(theta2)));

% Display the symbolic expression for the input force
disp('Symbolic expression for input force F_in:');
disp(F_in);

%% Assign Real Values to Parameters
% These are the physical parameters of the system
E = 190e9; % Young's modulus of Steel in Pascals (Pa)
l = 0.05; % Length of the compliant part (meters)
w = 0.05; % Width of the compliant part (meters)
t = 0.005; % Thickness of the compliant part (meters)

% Real values for lengths and spring constants
L2_real = 0.2; % Length L2 (meters)
L3_real = 0.25; % Length L3 (meters)
K1_real = (E * w * t^3) / (12 * l); % Spring constant K1 (N/rad)
K2_real = K1_real; % Assume K2 = K1
K3_real = K1_real; % Assume K3 = K1

% Initial angle values for theta2 and theta3
theta2_0_real = pi / 10; % Initial value of theta2 (radians)
theta3_0_real = 2 * pi - pi / 10; % Initial value of theta3 (radians)

% Substitute real values into the symbolic expression for F_in
F_in_real = subs(F_in, [K1, K2, K3, L2, L3, theta2_0, theta3_0], ...
                    [K1_real, K2_real, K3_real, L2_real, L3_real, theta2_0_real, theta3_0_real]);

%% Force Calculation for a Range of Angles
% Define a range of theta2 values for plotting the force
n = 101; % Number of points
theta2_real = linspace(theta2_0_real, pi/2, n); % Angle range for theta2 (in radians)
F_in_plot = zeros(n, 1); % Preallocate array for input forces

% Compute the input force for each value of theta2
for i = 1:length(theta2_real)
    F_in_plot(i) = double(subs(F_in_real, theta2, theta2_real(i))); % Calculate F_in at each theta2
end

% Convert theta2 to degrees for plotting
theta2_real_deg = theta2_real * 180 / pi;

% Convert force from Newtons to kilonewtons (kN)
F_in_plot_kN = F_in_plot / 1000;

%% Plotting the Input Force
% Plot the input force as a function of theta2
figure;
plot(theta2_real_deg, F_in_plot_kN, 'LineWidth', 2);
grid on;
xlabel('\theta_{2} [°]'); % Label for x-axis
ylabel('Input force required [kN]'); % Label for y-axis
xlim([0 100]); % Set x-axis limits
title('Input Force as a Function of \theta_2'); % Title for the plot

%% Potential Energy Calculation
% The potential energy is calculated by numerically integrating the input force
% over the range of theta2 values. The result is in kilojoules (kJ).
W = trapz(theta2_real_deg, F_in_plot_kN / 180 * pi) * L2_real * 2; % Energy in joules
W_kJ = W / 1000; % Convert to kilojoules

% Display the calculated potential energy
disp(['The potential energy of the jumping robot is: ', num2str(W_kJ), ' kJ']);

%% Stress Calculations
% The stress of the flextures are calculated numerically by using the
% Lagrange coordinates. The maximal stress occurs at the second lagrange
% coordinate. 

% Lagrange coordinate
phi_2 = (theta3 - theta3_0_real) - (theta2_real - theta2_0_real)

sigma_max = (E * t * phi_2) / (2 * L2);
