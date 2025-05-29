# Coordinates of the vertices
x1, y1 = 5, -2
x2, y2 = 10, 5
x3, y3 = 5, 5

# Calculate the area using the formula
area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

# Print the area rounded to the nearest tenth
print(round(area, 1))