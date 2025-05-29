# Coordinates of the points
x1, y1 = 0, 3  # D
x2, y2 = 0, 0  # B
x3, y3 = 8, 0  # C

# Calculate the area using the formula
area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
print(area)