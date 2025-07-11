import math

# The problem is to find the area of the smallest square enclosing a parametric curve.
# We first need to find the extent of the curve in x and y directions.

# x = cos(t) for t in [0, 2*pi). The range of cos(t) is [-1, 1].
x_min = -1
x_max = 1

# y is defined in two parts.
# For t in [0, pi], y = 1 * sin(t). The range of sin(t) is [0, 1].
y_max_part1 = 1
# For t in (pi, 2*pi), y = 5 * sin(t). The range of sin(t) is (-1, 0].
# So the range of 5*sin(t) is (-5, 0].
y_min_part2 = -5

# The total range for y is the union of the ranges from the two parts.
y_min = y_min_part2
y_max = y_max_part1

# Calculate the width and height of the bounding rectangle for the figure.
width = x_max - x_min
height = y_max - y_min

print(f"The figure's maximum width is {x_max} - ({x_min}) = {width}.")
print(f"The figure's maximum height is {y_max} - ({y_min}) = {height}.")

# The side of the smallest enclosing square is the maximum of the width and height.
side_of_square = max(width, height)
print(f"The side length of the smallest enclosing square is the larger dimension: {side_of_square}.")

# Calculate the area of the square.
area = side_of_square ** 2
print("The area of the square is found by squaring its side length.")
print(f"Final calculation: {int(side_of_square)} * {int(side_of_square)} = {int(area)}")
