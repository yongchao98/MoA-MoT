import math

# Step 1: Define the min and max for the x-coordinate from x = cos(t)
x_min = -1.0
x_max = 1.0

# Step 2: Define the min and max for the y-coordinate
# For t in [0, pi], y = sin(t), so y is in [0, 1]
# For t in (pi, 2pi), y = 5*sin(t), so y is in [-5, 0)
y_min = -5.0
y_max = 1.0

# Step 3: Calculate the width and height of the figure's bounding box
width = x_max - x_min
height = y_max - y_min

# Step 4: The side of the smallest outcircling square is the maximum of the width and height
side_length = max(width, height)

# Step 5: Calculate the area of the square
area = side_length ** 2

print(f"Width of the figure: {x_max} - ({x_min}) = {width}")
print(f"Height of the figure: {y_max} - ({y_min}) = {height}")
print(f"Side length of the smallest enclosing square: max({width}, {height}) = {side_length}")
print(f"The final equation for the area is:")
print(f"{int(side_length)} * {int(side_length)} = {int(area)}")
