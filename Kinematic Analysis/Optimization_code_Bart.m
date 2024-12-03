% Constants
clear all; clc;
%given in exam sheets:
gamma = 0.85;
K_theta = 2.65;
E = 190e9; % Young's modulus [Pa]           https://nl.china-stainless-steels.com/stainless-steel-plate/stainless-steel-1-4310.html
sigma_yield = 500e6; % Yield stress [Pa]    %given in lecture
theta20 = pi/2; % Initial angle in radians  
m = 0.1; % Mass [kg]                        % for now assumed but will follow from solidworks
g = 9.81; % Gravity [m/s^2]
h = 0.3; % Maximum height [m]                   % initial condition
x_final = 0.2; % Final horizontal position [m]  % initial condition

% Compute initial velocities and kinetic energy        % parabolic equation
% what initial velocity do we need to make the parabolic jump and how much energy does that generate??
v0y = sqrt(2 * g * h);                                      
t_up = v0y / g; % time it takes to move to the top of the parabolic
delta_h = h - 0.2;
t_down = sqrt(2 * delta_h / g); % time from top to 2/3 of the parabolic's height
t_total = t_up + t_down;
v0x = x_final / t_total;
v0 = sqrt(v0x^2 + v0y^2);
KE_initial = 0.5 * m * v0^2; % Total kinetic energy required [J]

% Adjust the required kinetic energy per mechanism
KE_per_mechanism = KE_initial / 2; % Energy per mechanism [J]

% Allowable thickness values in meters       % available values at IWS 
t_values = [0.0001, 0.00015, 0.0003, 0.0004];

% Step sizes            $ it needs to be fabricated so we want to use
% rounded values
step_w = 0.001; % 1 mm
step_L = 0.001; % 1 mm

% Variable ranges for w
w_min = 0.01; % Lower bound for w in meters
w_max = 0.05; % Upper bound for w in meters
w_values = w_min:step_w:w_max;

% Variable bounds for L
L_min = 0.02;
L_max = 0.06;
L_values = L_min:step_L:L_max;

% angle needs to be small enough so it does
% Variable ranges for delta_theta
delta_theta_min = deg2rad(1);
delta_theta_max = deg2rad(60);
step_theta = 0.01; % 0.01 rad

% Initialize variables to store optimal results
F_min_total = Inf;
best_energy_diff = Inf; % Initialize the best energy difference

% Loop over allowable thickness values
for t_idx = 1:length(t_values)
    t_fixed = t_values(t_idx);
    
    % Loop over all combinations of w and L
    for w = w_values
        for L = L_values
            % Calculate moment of inertia
            I = (1/12) * w * t_fixed^3;
            
            % Solve for delta_theta using the energy equality constraint
            delta_theta = sqrt( (KE_per_mechanism * L) / (4 * gamma * K_theta * E * I) );
            
            % Check if delta_theta is within bounds
            if delta_theta < delta_theta_min || delta_theta > delta_theta_max
                continue; % Skip this combination
            end
            
            % Discretize delta_theta to nearest step
            n_theta = round(delta_theta / step_theta);
            delta_theta_disc = n_theta * step_theta;
            
            % Check if discretized delta_theta is within bounds
            if delta_theta_disc < delta_theta_min || delta_theta_disc > delta_theta_max
                continue; % Skip this combination
            end

            % Recompute U with discretized delta_theta
            U = (4 * gamma * K_theta * E * I / L) * delta_theta_disc^2;
            energy_diff = abs(U - KE_per_mechanism);

            % Set a reasonable energy difference tolerance (e.g., 1% of KE_per_mechanism)
            energy_tolerance = 0.01 * KE_per_mechanism; % 1% tolerance
            if energy_diff > energy_tolerance
                continue; % Skip this combination
            end

            % Calculate maximum stress in the flexure
            sigma_max = (E * t_fixed * delta_theta_disc) / (2 * L);

            % Check stress constraint
            if sigma_max > sigma_yield
                continue; % Skip this combination
            end

            % Compute the force required per mechanism
            F = (4 * E * I * K_theta * delta_theta_disc) / (L^2 * cos(delta_theta_disc));

            % Check if this is the minimal force found so far
            if F < F_min_total
                % Store optimal variables
                F_min_total = F;
                t_opt = t_fixed;
                w_opt = w;
                L_opt = L;
                delta_theta_opt = delta_theta_disc;
                sigma_max_opt = sigma_max;
                U_opt = U;
                best_energy_diff = energy_diff;
            end
        end
    end
end

% Display optimal results
if isfinite(F_min_total)
    fprintf('Optimal Thickness t: %.6f m\n', t_opt);
    fprintf('Optimized Width w: %.6f m\n', w_opt);
    fprintf('Optimized Length L: %.6f m\n', L_opt);
    fprintf('Optimized Angular Deflection (theta2 - theta20): %.6f rad\n', delta_theta_opt);
    fprintf('Minimum Force Required per Mechanism: %.6f N\n', F_min_total);
    fprintf('Total Force Required: %.6f N\n', 2 * F_min_total);

    % Display constraints
    fprintf('Elastic Potential Energy U per Mechanism: %.6f J\n', U_opt);
    fprintf('Required Kinetic Energy per Mechanism: %.6f J\n', KE_per_mechanism);
    fprintf('Energy Difference: %.6f J\n', best_energy_diff);
    fprintf('Maximum Stress sigma_max: %.2f MPa\n', sigma_max_opt / 1e6);
else
    fprintf('No feasible solution found with the given parameters.\n');
end
