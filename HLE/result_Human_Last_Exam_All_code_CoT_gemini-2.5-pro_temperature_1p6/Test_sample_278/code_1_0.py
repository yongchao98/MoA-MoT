import math

# Step 1: Define the boundaries for the x and y coordinates based on the parametric equations.

# For x = cos(t), as t goes from 0 to 2*pi, the range of cos(t) is [-1, 1].
x_min = -1.0
x_max = 1.0

# For y, we analyze two parts:
# For t in [0, pi], y = 1 * sin(t). The maximum value is sin(pi/2) = 1.
y_max_part1 = 1.0
# For t in (pi, 2*pi), y = 5 * sin(t). The minimum value is 5*sin(3*pi/2) = -5.
y_min_part2 = -5.0

# The overall y range is the combination of the two parts.
y_min = y_min_part2
y_max = y_max_part1

# Step 2: Calculate the width and height of the bounding box of the figure.
width = x_max - x_min
height = y_max - y_min

# Step 3: The side of the smallest outcircling square is the maximum of the width and height.
side_length = max(width, height)

# Step 4: Calculate the area of the square.
area = side_length**2

# Print the results, including the final equation.
print(f"The full range of the x-coordinate is [{x_min}, {x_max}].")
print(f"The full range of the y-coordinate is [{y_min}, {y_max}].")
print(f"The width of the figure is {width}.")
print(f"The height of the figure is {height}.")
print(f"The side length of the smallest enclosing square is max({width}, {height}) = {side_length}.")
print("The area is calculated from the side length.")
print(f"The final equation is Area = {int(side_length)} * {int(side_length)}.")
print(f"The calculated area is {int(area)}.")
