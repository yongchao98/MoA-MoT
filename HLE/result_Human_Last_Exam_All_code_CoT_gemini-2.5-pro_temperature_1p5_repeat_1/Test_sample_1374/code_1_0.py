import math

# This script calculates the furthest distance from point A for a shape
# that maximizes the gravitational field at A.
# The formula for this distance is the cube root of (15 / (4 * pi)).

# Define the constants in the equation
numerator = 15
four_pi = 4 * math.pi

# Calculate the value inside the cube root
ratio = numerator / four_pi

# Calculate the cube root to find the maximum distance
max_distance = ratio**(1/3)

# Print the final equation with the result
print(f"The final equation is: ({numerator} / (4 * {math.pi}))^(1/3)")
print(f"The furthest point on the surface is {max_distance} meters from point A.")