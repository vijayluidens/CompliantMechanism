import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import data
file_path = r"C:\Users\vijay\WB - Masters\Compliant Mechanisms\Test Results\20mm.csv"
file_path_2 = r"C:\Users\vijay\WB - Masters\Compliant Mechanisms\Kinematic Analysis\Data.csv"
file_path_3 = r"C:\Users\vijay\WB - Masters\Compliant Mechanisms\APDL Code\APDL_DATA_1_flexture.csv"

# Load the CSV file with pandas
data = pd.read_csv(file_path, delimiter=",", skiprows=1)
data_2 = pd.read_csv(file_path_2, delimiter=",")
data_3 = pd.read_csv(file_path_3, delimiter=',')

# Convert data to numpy arrays
displacement = data.iloc[:, 0].to_numpy()  # First column
force = data.iloc[:, 1].to_numpy()         # Second column
prbm_displacement = data_2.iloc[:, 0].to_numpy()
prbm_force = data_2.iloc[:, 1].to_numpy()
apdl_displacement = data_3.iloc[:, 1].to_numpy()
apdl_force = data_3.iloc[:, 3].to_numpy()


# Determine the peak force, displacement & index
peak_index = np.argmax(force)
peak_force = np.max(force)
peak_displacement = np.max(displacement)

# PRBM Values
prbm_peak_force = np.max(prbm_force)
prbm_peak_displacement = np.max(prbm_displacement)

#APDL Values
apdl_peak_force = np.max(apdl_force)
apdl_peak_displacement = np.max(apdl_displacement)

# Ensure matching lengths for descending slice
descending_displacement = displacement[peak_index:]
descending_force = force[peak_index:]

if len(descending_displacement) > len(descending_force):
    descending_displacement = descending_displacement[:len(descending_force)]
elif len(descending_force) > len(descending_displacement):
    descending_force = descending_force[:len(descending_displacement)]

# Plot ascending
plt.figure(figsize=(10, 6))
plt.plot(displacement[:peak_index+1], force[:peak_index+1], label='Test: Load Phase', linewidth=5)

# Plot descending
plt.plot(descending_displacement, descending_force, linestyle='dotted', linewidth=4, label='Test: Deload Phase')

# Plot data from the prbm file
plt.plot(data_2['Deflection_X_mm'], data_2['Force_N'], label='PRBM Results', linewidth=4, linestyle = 'dashed')

# Plot data from the APDL file
plt.plot(data_3['Displacement (mm)'], data_3['Reaction_Force_2_Flextures (N)'], label='APDL Results', linewidth=4)

# Customize plot
plt.title('Force vs Displacement', fontsize=14)
plt.xlabel('Displacement (mm)', fontsize=12)
plt.ylabel('Force (N)', fontsize=12)
plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()

# Adjust x-axis and y-axis limits to add some padding
plt.xlim(0, peak_displacement + 2.5)  # Extend x-axis slightly beyond the last point
plt.ylim(0, peak_force + 6)  # Extend y-axis slightly beyond the last point


# Plot dot at maximum force value
plt.scatter(peak_displacement, peak_force, color='red', zorder=5)
plt.text(np.round(peak_displacement)-1, peak_force+1, f'({np.round(peak_displacement)}, {np.round(peak_force)})', fontsize=12)

# Plot dot at maximum force value: prbm
plt.scatter(prbm_peak_displacement, prbm_peak_force, color='blue', zorder=5)
plt.text(np.round(prbm_peak_displacement)-1, prbm_peak_force -2.5, f'({np.round(prbm_peak_displacement)}, {np.round(prbm_peak_force)})', fontsize=12)

# Plot dot at maximum force value: prbm
plt.scatter(apdl_peak_displacement, apdl_peak_force, color='green', zorder=5)
plt.text(np.round(apdl_peak_displacement)-1, apdl_peak_force +1, f'({np.round(apdl_peak_displacement)}, {np.round(apdl_peak_force)})', fontsize=12)

#Show Plot
plt.show()