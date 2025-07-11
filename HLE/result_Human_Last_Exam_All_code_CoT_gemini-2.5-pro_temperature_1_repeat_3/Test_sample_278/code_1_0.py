import numpy as np

# For t in [0, 2*pi], x = cos(t) ranges from -1 to 1.
x_min = -1
x_max = 1

# For t in [0, pi], y = 1*sin(t) ranges from 0 to 1.
# For t in (pi, 2*pi), y = 5*sin(t) ranges from -5 to 0.
# The total range for y is [-5, 1].
y_min = -5
y_max = 1

# Calculate the width and height of the bounding box
width = x_max - x_min
height = y_max - y_min

# The side of the smallest outcircling square is the maximum of the width and height
side_length = max(width, height)

# Calculate the area of the square
area = side_length * side_length

print(f"The range of x is [{x_min}, {x_max}], so the width is {x_max} - ({x_min}) = {width}.")
print(f"The range of y is [{y_min}, {y_max}], so the height is {y_max} - ({y_min}) = {height}.")
print(f"The side of the smallest square is the maximum of the width and height: max({width}, {height}) = {side_length}.")
print("The area of the square is side * side.")
print(f"Final calculation: {int(side_length)} * {int(side_length)} = {int(area)}")
