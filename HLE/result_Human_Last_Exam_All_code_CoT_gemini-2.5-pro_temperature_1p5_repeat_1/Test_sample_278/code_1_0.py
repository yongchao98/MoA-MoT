import math

# Step 1: Define the boundaries for x and y
# For x = cos(t), the range over [0, 2*pi] is [-1, 1].
x_min = -1
x_max = 1

# For y, we have two pieces.
# For t in [0, pi], y = sin(t). The range is [0, 1].
# For t in (pi, 2*pi), y = 5*sin(t). The range is [-5, 0).
# The total range for y is the union of these two, which is [-5, 1].
y_min = -5
y_max = 1

# Step 2: Calculate the width and height of the figure's bounding box
width = x_max - x_min
height = y_max - y_min

# Step 3: The side of the smallest enclosing square is the maximum of the width and height.
side_of_square = max(width, height)

# Step 4: Calculate the area of the square.
area = side_of_square * side_of_square

# Step 5: Print the results, including the final equation as requested.
print(f"The figure's width is {x_max} - ({x_min}) = {width}.")
print(f"The figure's height is {y_max} - ({y_min}) = {height}.")
print(f"The side of the smallest enclosing square must be at least {side_of_square}.")
print(f"The area of the smallest square is {side_of_square} * {side_of_square} = {area}.")
