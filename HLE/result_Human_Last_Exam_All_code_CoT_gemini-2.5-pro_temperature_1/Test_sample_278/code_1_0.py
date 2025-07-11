import math

# Step 1: Determine the range of the x-coordinate.
# x = cos(t) for t in [0, 2*pi].
# The range of cos(t) is [-1, 1].
x_min = -1
x_max = 1
print(f"The minimum value of x is {x_min}.")
print(f"The maximum value of x is {x_max}.")

# Step 2: Determine the range of the y-coordinate.
# For t in [0, pi], y = sin(t), which ranges from 0 to 1.
# For t in (pi, 2*pi), y = 5*sin(t), which ranges from -5 to 0.
# The overall range for y is the union of [0, 1] and [-5, 0), which is [-5, 1].
y_min = -5
y_max = 1
print(f"The minimum value of y is {y_min}.")
print(f"The maximum value of y is {y_max}.")

# Step 3: Calculate the width and height of the bounding box.
width = x_max - x_min
height = y_max - y_min
print(f"The total width of the figure is {x_max} - ({x_min}) = {width}.")
print(f"The total height of the figure is {y_max} - ({y_min}) = {height}.")

# Step 4: Determine the side length of the smallest enclosing square.
# The side length is the maximum of the width and height.
side_length = max(width, height)
print(f"The side length of the smallest enclosing square is max({width}, {height}) = {side_length}.")

# Step 5: Calculate the area of the square.
area = side_length ** 2
print(f"The area of the smallest enclosing square is {side_length} * {side_length} = {area}.")
