import math

# Step 1 & 2: Determine the overall range of x and y
# For x = cos(t), the range over [0, 2*pi] is [-1, 1].
x_min = -1
x_max = 1

# For y, we have two parts:
# Part 1: y = sin(t) for t in [0, pi]. Range is [0, 1].
# Part 2: y = 5*sin(t) for t in (pi, 2*pi). Range is [-5, 0).
# The overall range for y is the union of these two ranges.
y_min = -5
y_max = 1

# Step 3: Calculate the dimensions of the bounding box
width = x_max - x_min
height = y_max - y_min

# Step 4: Find the side length of the smallest enclosing square
side_length = max(width, height)

# Step 5: Calculate the area of the square
area = side_length * side_length

# Print the final calculation and result
print(f"The side length of the smallest enclosing square is the maximum of the figure's width and height.")
print(f"Width = {x_max} - ({x_min}) = {width}")
print(f"Height = {y_max} - ({y_min}) = {height}")
print(f"Side of the square = max({width}, {height}) = {side_length}")
print(f"The area of the square is:")
print(f"{int(side_length)} * {int(side_length)} = {int(area)}")
