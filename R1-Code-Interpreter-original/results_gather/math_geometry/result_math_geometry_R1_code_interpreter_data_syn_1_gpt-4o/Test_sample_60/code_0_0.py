import math

# Given values
radius = 10
distance_between_centers = 10

# Calculate the half-length of the chord using Pythagorean theorem
y = distance_between_centers / 2
x_squared = radius**2 - y**2
x = math.sqrt(x_squared)

# Full length of the chord
chord_length = 2 * x

# Output the chord length in simplest radical form
print(f"{chord_length:.2f}")