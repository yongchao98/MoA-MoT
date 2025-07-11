import math

# The curve is given by x = cos(t), y = b*sin(t)
# where b=1 for t in [0, pi] and b=5 for t in (pi, 2*pi).

# 1. Find the range of x = cos(t) for t in [0, 2*pi]
x_max = 1
x_min = -1

# 2. Find the range of y
# For t in [0, pi], y = sin(t). The range is [0, 1].
y_max_part1 = 1
y_min_part1 = 0

# For t in (pi, 2*pi), y = 5*sin(t). The range is [-5, 0).
y_max_part2 = 0
y_min_part2 = -5

# The overall range for y is the union of the ranges from the two parts.
y_max = max(y_max_part1, y_max_part2)
y_min = min(y_min_part1, y_min_part2)

# 3. Calculate the width and height of the bounding box
width = x_max - x_min
height = y_max - y_min

# 4. Determine the side of the smallest enclosing square
side_of_square = max(width, height)

# 5. Calculate the area of the square
area = side_of_square * side_of_square

print("To find the area of the smallest outcircling square, we first find the dimensions of the figure.")
print(f"The width of the figure is max(x) - min(x) = {x_max} - ({x_min}) = {width}")
print(f"The height of the figure is max(y) - min(y) = {y_max} - ({y_min}) = {height}")
print("")
print("The side of the smallest square that encloses the figure must be equal to the larger dimension.")
print(f"Side of square = max({width}, {height}) = {side_of_square}")
print("")
print("The area of this square is side * side.")
print(f"Area = {side_of_square} * {side_of_square} = {area}")

<<<36>>>