import numpy as np

# This script demonstrates how the tanh function can be used to cluster grid points,
# a common technique for controlling grid resolution in computational modeling.

# 1. Define the parameters for the transformation equation
# A higher stretching_factor leads to more pronounced clustering around the center.
stretching_factor = 2.5
num_points = 11

# 2. Create a uniform grid in computational space (from -1 to 1)
uniform_points = np.linspace(-1.0, 1.0, num_points)

# 3. Apply the tanh transformation and print the details for each step.
# The equation is: new_point = tanh(stretching_factor * old_point) / tanh(stretching_factor)
# We will calculate and print each component of this equation.

print("Demonstrating grid clustering using the tanh function.")
print("The transformation equation is: new_point = tanh(A * x) / tanh(A)")
print(f"Where A (stretching_factor) = {stretching_factor}\n")
print("Calculation for each point 'x' in the uniform grid:")
print("-----------------------------------------------------------------")
print("  x (old) |   A*x    | tanh(A*x) |  tanh(A)  |  x (new)  |")
print("-----------------------------------------------------------------")

# Calculate the denominator once, as it's constant for this transformation
denominator = np.tanh(stretching_factor)

# Loop through each point in the uniform grid
for x in uniform_points:
    # Calculate each part of the equation
    a_times_x = stretching_factor * x
    numerator = np.tanh(a_times_x)
    new_x = numerator / denominator

    # Print the numbers used in the equation for this specific point
    print(f"{x:9.4f} | {a_times_x:8.4f} | {numerator:9.4f} | {denominator:9.4f} | {new_x:9.4f} |")
print("-----------------------------------------------------------------")
print("\nNotice how the 'x (new)' points are clustered more densely around 0.0.")
