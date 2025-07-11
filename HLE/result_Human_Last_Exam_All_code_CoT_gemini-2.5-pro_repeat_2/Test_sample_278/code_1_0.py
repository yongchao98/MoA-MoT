import math

# Step 1: Find the minimum and maximum x and y values.

# For x = cos(t), the range over [0, 2*pi) is [-1, 1].
x_min = -1.0
x_max = 1.0

# For y, we have two parts.
# Part 1: y = 1 * sin(t) for t in [0, pi]. Range is [0, 1].
# Part 2: y = 5 * sin(t) for t in (pi, 2*pi). Range is [-5, 0).
# The overall range for y is the union of these two ranges.
y_min = -5.0
y_max = 1.0

# Step 2: Calculate the width and height of the bounding box.
width = x_max - x_min
height = y_max - y_min

# Step 3: The side of the smallest enclosing square is the maximum of the width and height.
side_length = max(width, height)

# Step 4: Calculate the area of the square.
area = side_length ** 2

# Print the calculation steps
print(f"The figure spans from x={x_min} to x={x_max}, so its width is {x_max} - ({x_min}) = {width}.")
print(f"The figure spans from y={y_min} to y={y_max}, so its height is {y_max} - ({y_min}) = {height}.")
print(f"The side of the smallest enclosing square must be the larger of the width and height, which is max({width}, {height}) = {side_length}.")
print(f"The area of this square is {side_length} * {side_length} = {area}.")