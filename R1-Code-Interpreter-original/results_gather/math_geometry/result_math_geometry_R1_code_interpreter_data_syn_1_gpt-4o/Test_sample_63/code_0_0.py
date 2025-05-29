import math

# Given values
theta = 252  # degrees
r = 10  # radius of the circle

# Calculate the arc length
arc_length = (theta / 360) * 2 * math.pi * r

# Calculate the base radius of the cone
base_radius = arc_length / (2 * math.pi)

# Print the base radius
print(base_radius)