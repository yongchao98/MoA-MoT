import numpy as np

# Step 1: Define the ranges for x and y based on the parametric equations.

# For x = cos(t), the range over t in [0, 2*pi) is [-1, 1].
x_min = -1
x_max = 1

# For y = b*sin(t):
# When t is in [0, pi], b=1, so y = sin(t) ranges from 0 to 1.
# When t is in (pi, 2*pi), b=5, so y = 5*sin(t) ranges from -5 to 0.
# The combined range for y is therefore [-5, 1].
y_min = -5
y_max = 1

# Step 2: Calculate the dimensions of the bounding box.
width = x_max - x_min
height = y_max - y_min

print(f"The overall range for x is [{x_min}, {x_max}].")
print(f"The overall range for y is [{y_min}, {y_max}].")
print(f"The width of the bounding box for the figure is {x_max} - ({x_min}) = {width}.")
print(f"The height of the bounding box for the figure is {y_max} - ({y_min}) = {height}.")

# Step 3: Determine the side length of the smallest enclosing square.
# The side length must be the maximum of the width and height.
side_length = max(width, height)
print(f"\nThe side of the smallest square that outcircles the figure is the maximum of its width and height.")
print(f"Side length = max({width}, {height}) = {side_length}.")

# Step 4: Calculate the area of the square.
area = side_length ** 2
print(f"\nThe area of the square is the side length squared.")
print(f"Area = {side_length} * {side_length} = {area}")
