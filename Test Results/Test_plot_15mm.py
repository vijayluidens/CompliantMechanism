import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read data from a CSV file, skipping the first row
file_path = r"C:\Users\vijay\WB - Masters\Compliant Mechanisms\Test Results\15mm_test.csv"  

# Load the CSV file with pandas
data = pd.read_csv(file_path, delimiter=",", skiprows=1)

# Convert data to numpy arrays
displacement = data.iloc[:, 0].to_numpy()  # First column
force = data.iloc[:, 1].to_numpy()         # Second column

# Determine the peak force index
peak_index = np.argmax(force)

# Ensure matching lengths for descending slice
descending_displacement = displacement[peak_index:]
descending_force = force[peak_index:]

if len(descending_displacement) > len(descending_force):
    descending_displacement = descending_displacement[:len(descending_force)]
elif len(descending_force) > len(descending_displacement):
    descending_force = descending_force[:len(descending_displacement)]

# Plot ascending
plt.figure(figsize=(10, 6))
plt.plot(displacement[:peak_index+1], force[:peak_index+1], label='Ascending', linewidth=2)

# Plot descending
plt.plot(descending_displacement, descending_force, linestyle='dotted', linewidth=2, label='Descending')

# Customize plot
plt.title('Force vs Displacement', fontsize=14)
plt.xlabel('Displacement (mm)', fontsize=12)
plt.ylabel('Force (N)', fontsize=12)
plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()

# Show plot
plt.show()
