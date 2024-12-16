import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp

# Import data
file_path = r"C:\Users\vijay\WB - Masters\Compliant Mechanisms\Test Results\20mm.csv"

# Load the CSV file with pandas
data = pd.read_csv(file_path, delimiter=",", skiprows=1)

# Convert data to numpy arrays
displacement = data.iloc[:, 0].to_numpy()  # First column
force = data.iloc[:, 1].to_numpy()         # Second column

# Determine the peak force, displacement & index
peak_index = np.argmax(force)
peak_force = np.max(force)
peak_displacement = np.max(displacement)

# Ensure matching lengths for descending slice
descending_displacement = displacement[peak_index:]
descending_force = force[peak_index:]

if len(descending_displacement) > len(descending_force):
    descending_displacement = descending_displacement[:len(descending_force)]
elif len(descending_force) > len(descending_displacement):
    descending_force = descending_force[:len(descending_displacement)]

# Plot ascending
plt.figure(figsize=(10, 6))
plt.plot(displacement[:peak_index+1], force[:peak_index+1], label='Ascending', linewidth=5)

# Plot descending
plt.plot(descending_displacement, descending_force, linestyle='dotted', linewidth=4, label='Descending')

# Customize plot
plt.title('Force vs Displacement', fontsize=14)
plt.xlabel('Displacement (mm)', fontsize=12)
plt.ylabel('Force (N)', fontsize=12)
plt.grid(True)
plt.legend(fontsize=12)
# plt.tight_layout()

# Plot dot at maximum force value
plt.scatter(peak_displacement, peak_force, color='red', zorder=5)
plt.text(np.round(peak_displacement), peak_force, f'({np.round(peak_displacement)}, {np.round(peak_force):.4f})', fontsize=12, verticalalignment='top')


# Show plot
plt.show()

# Define symbolic variables
X = sp.symbols('X', real=True)
K1, K2, K3 = sp.symbols('K1 K2 K3', real=True)

# Define parameters
L = 125  # Length of the link in mm
L_flexure = 50  # Length of the flexure
E = 200  # Young's modulus

t = 0.4  # Thickness in mm
w = 30  # Width in mm
theta_20 = np.deg2rad(18)  # Initial angle theta_2 in radians (18 degrees)
theta_30 = 2 * np.pi - theta_20  # Initial angle theta_3
X_initial = 2 * L * np.cos(theta_20)  # Initial value of X based on theta_20

# Define the range of X within valid bounds
X_range = np.linspace(X_initial, X_initial - 20, 10)

# Calculate theta_2 and theta_3 as functions of X
theta_2_formula = sp.acos(X / (2 * L))
theta_3_formula = 2 * sp.pi - theta_2_formula

# Calculate d(theta_2)/dX and d(theta_3)/dX:
dtheta_2_dX = sp.diff(theta_2_formula, X)
dtheta_3_dX = sp.diff(theta_3_formula, X)

# Define phi values relative to initial conditions
phi1 = theta_2_formula - theta_20
phi2 = (theta_3_formula - theta_30) - (theta_2_formula - theta_20)
phi3 = theta_3_formula - theta_30

# Calculate spring constants
K = (E) * (t * w**3 / 12) / (L_flexure)  # Simplified flexure stiffness equation
K1 = K2 = K3 = K

# Force expression
Force_expr = (
    phi1 * K1 * dtheta_2_dX +
    phi2 * K2 * (dtheta_3_dX - dtheta_2_dX) +
    phi3 * K3 * dtheta_3_dX
)

# Calculate the initial force when X = X_initial
Initial_force = float(sp.N(Force_expr.subs(X, X_initial)))

# Pre-allocate arrays
theta_2_values = []
theta_3_values = []
force_values = []

# Loop through each X value
for X_value in X_range:
    Force_sub_X = Force_expr.subs(X, X_value)

    # Evaluate theta_2 and theta_3
    theta_2_values.append(float(sp.N(theta_2_formula.subs(X, X_value))))
    theta_3_values.append(float(sp.N(theta_3_formula.subs(X, X_value))))

    # Evaluate force and subtract the initial force
    force_values.append(float(sp.N(Force_sub_X)) - Initial_force)

# Convert results to numpy arrays for plotting
deflection_X = -(X_range - X_initial)  # Inverted displacement values
force_values = np.array(force_values)
theta_2_values = np.array(theta_2_values)
theta_3_values = np.array(theta_3_values)

# Plot Force-Deflection graph
plt.figure()
plt.plot(deflection_X, force_values, '--sm')  # Invert the displacement values
plt.xlabel('Deflection X (mm)')
plt.ylabel('Force (N)')
plt.title('Force-Deflection Graph')
plt.ylim([0, None])  # Set y-axis minimum to 0
plt.grid()
plt.show()

# Display results
print('Theta_2 values (deg):')
print(np.rad2deg(theta_2_values))
print('Theta_3 values (deg):')
print(np.rad2deg(theta_3_values))
print('Force values (N):')
print(force_values)

# Export data to an Excel file
plot_data = pd.DataFrame({
    'Deflection_X_mm': deflection_X,
    'Force_N': force_values,
    'Theta2_deg': np.rad2deg(theta_2_values),
    'Theta3_deg': np.rad2deg(theta_3_values)
})
filename = 'Force_Deflection_Data.xlsx'
plot_data.to_excel(filename, index=False)
