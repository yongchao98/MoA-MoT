import math

# Step 1: Determine the range of the x-coordinates.
# x = cos(t) for t in [0, 2pi).
x_min = -1
x_max = 1
width = x_max - x_min

# Step 2: Determine the range of the y-coordinates.
# For t in [0, pi], y = sin(t), so y is in [0, 1].
# For t in (pi, 2pi), y = 5*sin(t), so y is in [-5, 0).
# The total range for y is the union of these two intervals.
y_min = -5
y_max = 1
height = y_max - y_min

# Step 3: The side of the smallest enclosing square is the maximum of the width and height.
side_length = max(width, height)

# Step 4: Calculate the area of the square.
area = side_length * side_length

# Step 5: Print the final equation for the area.
print(f"The width of the figure is {width}.")
print(f"The height of the figure is {height}.")
print(f"The side length of the smallest enclosing square is {side_length}.")
print(f"The area of the square is {int(side_length)} * {int(side_length)} = {int(area)}")
