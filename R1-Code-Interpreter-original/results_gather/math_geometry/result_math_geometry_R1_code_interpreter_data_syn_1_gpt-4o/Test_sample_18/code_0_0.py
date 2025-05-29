import math

# Given values
radius = 4
diameter = 2 * radius

# Calculate the side of the square
side_of_square = diameter / math.sqrt(2)

# Calculate the area of the square
area_of_square = side_of_square ** 2

# Calculate the area of the circle
area_of_circle = math.pi * radius ** 2

# Calculate the area of the shaded region
shaded_area = area_of_circle - area_of_square

print(shaded_area)