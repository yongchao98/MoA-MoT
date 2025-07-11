import numpy as np

# This script demonstrates how the hyperbolic tangent (tanh) function
# is used for grid generation to control spacing and ensure resolution consistency.
# It works by mapping a uniform grid to a non-uniform, clustered grid.

# The general mapping equation is: x_physical = L * tanh(c * x_computational) / tanh(c)
# L: The half-width of the physical domain.
# c: A clustering parameter. c > 1 creates clustering around the center.
# x_computational: The original, uniformly spaced coordinates, typically from -1 to 1.
# x_physical: The resulting, non-uniformly spaced coordinates.

# --- Parameters for our demonstration ---
L = 1.0
c = 2.5
num_points = 11

# The final equation we will demonstrate is:
# x_physical = 1.0 * tanh(2.5 * x_computational) / tanh(2.5)
# Here, we output the numbers used in this specific equation.
print("Demonstration of grid clustering using the tanh function.")
print("The mapping equation is: x_physical = L * tanh(c * x_computational) / tanh(c)")
print(f"The specific equation used here is: x_physical = {L} * tanh({c} * x_computational) / tanh({c})")
print("-" * 50)

# 1. Create a uniform grid in the computational domain.
x_computational = np.linspace(-1.0, 1.0, num_points)
print("Original uniform grid points:")
print(np.round(x_computational, 4))
print("-" * 50)

# 2. Apply the tanh mapping function to create the clustered physical grid.
denominator = np.tanh(c)
x_physical = L * np.tanh(c * x_computational) / denominator
print("New grid points clustered around the center:")
print(np.round(x_physical, 4))
print("-" * 50)
print("As you can see, the points in the new grid are more densely packed near 0.0.")
print("This is how tanh is used to increase resolution in a specific region.")
