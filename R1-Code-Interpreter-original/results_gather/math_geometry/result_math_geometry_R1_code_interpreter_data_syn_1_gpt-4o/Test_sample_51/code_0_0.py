import sympy as sp

# Define the height of the equilateral triangle
side_length = 8
height = (sp.sqrt(3) / 2) * side_length

# Calculate the y-coordinate of the third vertex
y_coordinate = 5 + height

# Simplify the expression
simplified_y_coordinate = sp.simplify(y_coordinate)

print(simplified_y_coordinate)