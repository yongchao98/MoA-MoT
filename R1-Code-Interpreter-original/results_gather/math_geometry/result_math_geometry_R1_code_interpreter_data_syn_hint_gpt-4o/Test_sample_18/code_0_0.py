import math

# Given values
radius = 4
side_of_square = 4 * math.sqrt(2)

# Calculate areas
area_of_circle = math.pi * radius**2
area_of_square = side_of_square**2

# Calculate shaded area
shaded_area = area_of_circle - area_of_square

print(shaded_area)