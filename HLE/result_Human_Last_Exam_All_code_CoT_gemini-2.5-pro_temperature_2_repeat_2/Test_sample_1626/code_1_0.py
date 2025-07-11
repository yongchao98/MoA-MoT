import numpy as np

# This script demonstrates how the hyperbolic tangent (tanh) function
# is used to create non-uniform grid spacing from a uniform set of points.
# This technique is essential for achieving resolution consistency in
# complex geometries like a torus.

# 1. Define a set of uniformly spaced points in a computational coordinate system.
# Let's create 11 points from -3 to 3.
uniform_points = np.linspace(-3.0, 3.0, 11)
print("Original uniform grid points:")
print(np.round(uniform_points, 4))
print("-" * 50)

# 2. Apply the tanh function to transform these points.
# The transformation equation is: new_coord = tanh(old_coord)
# The tanh function maps the points to the range (-1, 1), clustering them
# around 0.
transformed_points = np.tanh(uniform_points)
print("Transformed grid points using tanh:")
print("Notice how the points are now clustered around 0.0")
print(np.round(transformed_points, 4))
print("-" * 50)

# 3. Fulfilling the requirement to output the equation with each number.
# We will show the equation for the first three points as an example.
print("The final equation is of the form: transformed_point = tanh(original_point)")
print("Showing the equation with numbers for the first few points:")
for i in range(len(uniform_points)):
    original_val = uniform_points[i]
    transformed_val = transformed_points[i]
    # Here we output each number (the result and the input) in the final equation
    print(f"{transformed_val:.4f} = tanh({original_val:.4f})")
