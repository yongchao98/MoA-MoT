import math

# Given values
radius = 3
height = 10

# Calculate lateral surface area
lateral_surface_area = 2 * math.pi * radius * height

# Calculate area of the two ends
ends_area = 2 * math.pi * radius**2

# Calculate total surface area
total_surface_area = lateral_surface_area + ends_area

# Output the total surface area
print(total_surface_area)