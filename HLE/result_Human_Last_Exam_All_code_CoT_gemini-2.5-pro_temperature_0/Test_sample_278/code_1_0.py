import math

# 1. Define the ranges for x and y
# For x = cos(t), as t goes from 0 to 2*pi, the range is [-1, 1]
x_min = -1
x_max = 1

# For y, we have two cases:
# Case 1: t in [0, pi], y = sin(t). Range is [0, 1].
# Case 2: t in (pi, 2*pi), y = 5*sin(t). Range is [-5, 0).
# The combined range for y is [-5, 1].
y_min = -5
y_max = 1

# 2. Calculate the width and height of the bounding box
width = x_max - x_min
height = y_max - y_min

# 3. The side of the smallest enclosing square is the maximum of the width and height
side_length = max(width, height)

# 4. Calculate the area of the square
area = side_length ** 2

print(f"The width of the figure is {x_max} - ({x_min}) = {width}")
print(f"The height of the figure is {y_max} - ({y_min}) = {height}")
print(f"The side length of the smallest enclosing square is max({width}, {height}) = {side_length}")
print(f"The area of the square is {side_length} * {side_length} = {area}")
<<<36>>>