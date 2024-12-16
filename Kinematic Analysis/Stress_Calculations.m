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

% Evaluate stresses over the range of X
sigma_1_values = zeros(size(X_range));
sigma_2_values = zeros(size(X_range));
sigma_3_values = zeros(size(X_range));

for i = 1:length(X_range)
    X_value = X_range(i);

    % Substitute X into the updated phi relations
    phi1_current = subs(phi1, X, X_value);
    phi2_current = subs(phi2, X, X_value);
    phi3_current = subs(phi3, X, X_value);

    % Update moments
    M_1_current = double(subs(M_1, X, X_value));
    M_2_current = double(subs(M_2, X, X_value));
    M_3_current = double(subs(M_3, X, X_value));

    % Calculate stresses
    sigma_1_values(i) = double(M_1_current * C / I_value);
    sigma_2_values(i) = double(M_2_current * C / I_value);
    sigma_3_values(i) = double(M_3_current * C / I_value);
end

% Plot stress vs. Theta_2
figure;
plot(rad2deg(theta_2_values), abs(sigma_1_values), 'r', 'DisplayName', '\sigma_1');
hold on;
plot(rad2deg(theta_2_values), sigma_2_values, 'g', 'DisplayName', '\sigma_2');
plot(rad2deg(theta_2_values), sigma_3_values, 'b', 'DisplayName', '\sigma_3');
xlabel('Theta_2 (degrees)');
ylabel('Stress (Pa)');
title('Stress vs. Theta_2');
legend('show');
grid on;