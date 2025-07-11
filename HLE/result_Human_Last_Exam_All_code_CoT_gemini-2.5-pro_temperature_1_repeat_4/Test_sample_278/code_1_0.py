import math

# Define the boundaries for x and y
# For t in [0, pi], x=cos(t) and y=sin(t).
# x ranges from cos(pi)=-1 to cos(0)=1.
# y ranges from sin(0)=0 to sin(pi/2)=1.
x_max_1 = 1
x_min_1 = -1
y_max_1 = 1
y_min_1 = 0

# For t in (pi, 2pi), x=cos(t) and y=5*sin(t).
# x ranges from cos(pi)=-1 to cos(2pi)=1.
# y ranges from 5*sin(3pi/2)=-5 to a value approaching 5*sin(2pi)=0.
x_max_2 = 1
x_min_2 = -1
y_max_2 = 0
y_min_2 = -5

# Find the overall min and max values for the entire figure
x_max = max(x_max_1, x_max_2)
x_min = min(x_min_1, x_min_2)
y_max = max(y_max_1, y_max_2)
y_min = min(y_min_1, y_min_2)

# Calculate the width and height of the bounding box
width = x_max - x_min
height = y_max - y_min

# The side of the smallest enclosing square is the max of width and height
side_length = max(width, height)

# Calculate the area of the square
area = side_length ** 2

print(f"The bounding box of the figure has:")
print(f"Width = {x_max} - ({x_min}) = {width}")
print(f"Height = {y_max} - ({y_min}) = {height}")
print(f"\nThe side of the smallest enclosing square must be the larger of these two dimensions, which is {side_length}.")
print(f"The area of this square is {int(side_length)} * {int(side_length)} = {int(area)}.")
