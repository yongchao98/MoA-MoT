import numpy as np

# This script demonstrates how the tanh function is used to create
# non-uniform grid spacing to improve resolution consistency.

# Number of grid points to generate
num_points = 11

# --- Case 1: Uniform Spacing ---
# This is the naive approach, which leads to resolution issues on a torus.
# We create a set of points uniformly distributed between 0 and 360 degrees.
print("--- 1. Uniform Grid Spacing ---")
uniform_angles = np.linspace(0.0, 360.0, num_points)
print("A uniform grid results in these angles (in degrees):")
# We output each number in the resulting array.
print(np.round(uniform_angles, 2))
print("\n" + "="*50 + "\n")


# --- Case 2: Tanh-based Spacing ---
# This is the standard method for controlling grid packing.
# The mathematical function at the core of this transformation is tanh.
print("--- 2. Tanh-based Grid Spacing ---")
# 'packing_strength' controls how much the points are clustered at the ends.
# A value of 0 would be uniform; higher values increase clustering.
packing_strength = 2.5

# Create a symmetric input for the tanh function from -strength to +strength.
# The equation for this input is: x = linspace(-s, s, N)
tanh_input = np.linspace(-packing_strength, packing_strength, num_points)

# Apply the tanh function and then normalize the result to a [0, 1] range.
# The equation for the normalized spacing is: y = (tanh(x) + 1) / 2
tanh_normalized = (np.tanh(tanh_input) + 1) / 2.0

# Map the non-uniform [0, 1] spacing to the full 0-360 degree range.
tanh_angles = tanh_normalized * 360.0

print("Using the tanh function allows us to cluster grid points.")
print("The final equation for the angles is: angle = ((tanh(x) + 1) / 2) * 360")
print("This results in these non-uniformly spaced angles (in degrees):")
# We output each number in the final array.
print(np.round(tanh_angles, 2))
print("\nNotice how the points are now clustered near 0 and 360.")
