import math

# Given volume
volume = 36 * math.pi

# Calculate radius
r_cubed = (3 * 36) / 4
radius = r_cubed ** (1/3)

# Calculate surface area
surface_area = 4 * math.pi * radius**2

print(surface_area)