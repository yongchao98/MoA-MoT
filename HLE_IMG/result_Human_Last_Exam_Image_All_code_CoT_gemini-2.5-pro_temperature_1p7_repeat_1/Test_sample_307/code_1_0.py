import numpy as np

# Step 1: Define the data
# The infinity-norms for the nine Wigner 3-j symbols shown in the image,
# sourced from Table 1 of "Visualizing Wigner 3-j symbols" by M.J.D. Powell.
norms = [
    0.2222,  # Plot 1
    0.1601,  # Plot 2
    0.2795,  # Plot 3
    0.2222,  # Plot 4
    0.1601,  # Plot 5
    0.1068,  # Plot 6
    0.1601,  # Plot 7
    0.1068,  # Plot 8
    0.2795,  # Plot 9
]

# Step 2: Find the maximum and minimum values in the list
max_norm = max(norms)
min_norm = min(norms)

# Step 3: Calculate the ratio
ratio = max_norm / min_norm

# Step 4: Print the full calculation and result
print(f"The list of the nine infinity-norms is: {norms}")
print(f"Maximum infinity-norm = {max_norm}")
print(f"Minimum infinity-norm = {min_norm}")
print(f"Ratio = Maximum / Minimum")
print(f"Ratio = {max_norm} / {min_norm}")
print(f"The calculated ratio is: {ratio}")