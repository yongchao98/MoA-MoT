import numpy as np

# This script demonstrates how the hyperbolic tangent (tanh) function
# is used to create a non-uniform grid by "stretching" a uniform one.
# This technique is a cornerstone of controlling grid resolution in
# computational physics and engineering.

# 1. Define parameters for our 1D grid example.
num_points = 11  # The total number of points in our grid.
stretching_factor = 2.5  # A factor to control the intensity of the clustering. Higher means more clustering.

# 2. Create a uniformly spaced grid in what is called "computational space".
# These points are perfectly equidistant.
uniform_grid = np.linspace(-1.0, 1.0, num_points)

# 3. Use the tanh function to map the uniform grid to a new, non-uniform grid.
# The core of the transformation is the tanh function.
# We apply it and then rescale the output to fit back into the original [-1, 1] domain
# for a clear comparison.
stretched_grid_unscaled = np.tanh(stretching_factor * uniform_grid)
max_val = np.tanh(stretching_factor) # The value at the boundary
stretched_grid = stretched_grid_unscaled / max_val

# 4. Print the original and new grid point coordinates to show the effect.
# The printed numbers represent the coordinates of each point on a line.
print("Original Uniform Grid Points:")
# This loop outputs each number in the original set of coordinates.
for point in uniform_grid:
    print(f"{point:9.5f}")

print("\n" + "="*40 + "\n")

print(f"Stretched Grid Points (using tanh with factor={stretching_factor}):")
# This loop outputs each number in the final transformed set of coordinates.
for point in stretched_grid:
    print(f"{point:9.5f}")

print("\n" + "="*40 + "\n")
print("Notice how the points in the stretched grid are now clustered around 0.0,")
print("with the spacing increasing towards the ends at -1.0 and 1.0.")
