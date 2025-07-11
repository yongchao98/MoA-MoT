import math

# 1. Define the ranges for x and y
# For x = cos(t), t in [0, 2*pi), the range is [-1, 1].
x_min = -1.0
x_max = 1.0

# For y, we analyze two parts:
# Part 1: t in [0, pi], y = 1 * sin(t). Range of sin(t) is [0, 1].
y_max_1 = 1.0
y_min_1 = 0.0

# Part 2: t in (pi, 2*pi), y = 5 * sin(t). Range of sin(t) is [-1, 0).
# So range of 5*sin(t) is [-5, 0).
y_max_2 = 0.0
y_min_2 = -5.0

# Overall y range is the union of the two parts.
y_min = min(y_min_1, y_min_2)
y_max = max(y_max_1, y_max_2)

# 2. Calculate the width and height of the figure's bounding box.
width = x_max - x_min
height = y_max - y_min

# 3. Determine the side length of the smallest enclosing square.
# The side must be the maximum of the width and height.
side_length = max(width, height)

# 4. Calculate the area of the square.
area = side_length ** 2

# 5. Print the result.
# The final equation is side * side = area. We print each number.
print(f"The side length of the smallest enclosing square is {int(side_length)}.")
print(f"The area is calculated as: {int(side_length)} * {int(side_length)} = {int(area)}")
