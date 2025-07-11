import math

# Define the boundaries for the two parts of the curve
# Part 1: t in [0, pi], x=cos(t), y=sin(t)
x_min1 = -1.0
x_max1 = 1.0
y_min1 = 0.0
y_max1 = 1.0

# Part 2: t in (pi, 2pi), x=cos(t), y=5*sin(t)
x_min2 = -1.0
x_max2 = 1.0
y_min2 = -5.0
y_max2 = 0.0

# Determine the overall minimum and maximum coordinates for the entire figure
x_min_total = min(x_min1, x_min2)
x_max_total = max(x_max1, x_max2)
y_min_total = min(y_min1, y_min2)
y_max_total = max(y_max1, y_max2)

# Calculate the width and height of the figure's bounding box
width = x_max_total - x_min_total
height = y_max_total - y_min_total

# The side length of the smallest enclosing square is the maximum of the width and height
side_length = max(width, height)

# Calculate the area of the square
area = side_length ** 2

# Print the step-by-step calculation
print(f"The total width of the figure is {x_max_total} - ({x_min_total}) = {width}")
print(f"The total height of the figure is {y_max_total} - ({y_min_total}) = {height}")
print(f"The side of the smallest square that outcircles the figure must be the maximum of its width and height.")
print(f"Side length = max({width}, {height}) = {side_length}")
print(f"The area of the square is {side_length} * {side_length} = {area}")