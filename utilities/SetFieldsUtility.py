import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# Set number of levels within leaf area density field and starting and ending heights
Noincrements= 100
start_z = 0.5
end_z = 20
increment = (end_z - start_z) / Noincrements

# Set starting and end values for leaf area density
min_lad = 0.05
max_lad = 0.3

# Calculate center of the Gaussian
center_z = 13  # peak is at z = 13

# Objective function for the optimizer to minimize
# We want to find the width that gives LAD = 0.05 at start_z and LAD = 0.1 at end_z
def objective(width):
    start_lad = max_lad * np.exp(-0.5 * ((start_z - center_z) / width) ** 2)
    end_lad = max_lad * np.exp(-0.5 * ((end_z - center_z) / width) ** 2)
    return (start_lad - min_lad) ** 2 + (end_lad - 0.1) ** 2

# Use scipy's minimize function to find the optimal width
result = minimize(objective, (end_z - start_z) / 2)
width_z = result.x[0]

# Initialize current values
current_z = start_z

lad_values = []
z_values = [] 

with open('box_definitions.txt', 'w') as f:
    f.write("regions\n")
    f.write("(\n")

    for i in range(Noincrements):
        # Calculate LAD based on Gaussian function
        current_lad = max_lad * np.exp(-0.5 * ((current_z - center_z) / width_z) ** 2)

        # Store values for plot
        lad_values.append(current_lad)
        z_values.append(current_z)

        # Write boxToCell definition to file
        f.write("    boxToCell\n")
        f.write("    {\n")
        f.write(f"        box (200 0 {format(current_z, '.3f')}) (800 1 {format(current_z+increment, '.3f')});\n")
        f.write("        fieldValues\n")
        f.write("        (\n")
        f.write(f"            volScalarFieldValue LAD {format(current_lad, '.3f')} \n")
        f.write("            volScalarFieldValue Cd 0.2\n")
        f.write("        );\n")
        f.write("    }\n")
        
        # Update z position
        current_z += increment

    f.write(");")

# Generate plot
plt.plot(z_values, lad_values)
plt.xlabel('z position')
plt.ylabel('LAD')
plt.title('z vs LAD')
plt.grid(True)
plt.savefig("lad_plot.png")
plt.show()
