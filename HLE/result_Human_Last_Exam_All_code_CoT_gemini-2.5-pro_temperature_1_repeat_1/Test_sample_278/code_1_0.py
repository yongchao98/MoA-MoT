import math

# Step 1: Define the bounds of the figure based on the parametric equations.
# For x = cos(t), the minimum value is -1 and the maximum is 1 for t in [0, 2*pi].
x_min = -1.0
x_max = 1.0

# For y, we consider two parts:
# Part 1: t in [0, pi], y = sin(t). Max y is sin(pi/2) = 1. Min y is 0.
# Part 2: t in (pi, 2*pi), y = 5*sin(t). Min y is 5*sin(3*pi/2) = -5. Max y is 0.
# The overall y_min is the minimum of {0, -5}, and y_max is the maximum of {1, 0}.
y_min = -5.0
y_max = 1.0

print(f"The minimum and maximum x-values are {x_min} and {x_max}.")
print(f"The minimum and maximum y-values are {y_min} and {y_max}.")
print("-" * 30)

# Step 2: Calculate the width and height of the bounding box.
width = x_max - x_min
height = y_max - y_min

print(f"The width of the figure is {x_max} - ({x_min}) = {width}")
print(f"The height of the figure is {y_max} - ({y_min}) = {height}")
print("-" * 30)

# Step 3: The side of the smallest enclosing square is the maximum of the width and height.
side_length = max(width, height)

print(f"The side length of the smallest enclosing square is max({width}, {height}) = {side_length}")
print("-" * 30)

# Step 4: Calculate the area of the square.
area = side_length ** 2

print("The area of the square is side * side.")
print(f"Final Area = {side_length} * {side_length} = {area}")

<<<36.0>>>