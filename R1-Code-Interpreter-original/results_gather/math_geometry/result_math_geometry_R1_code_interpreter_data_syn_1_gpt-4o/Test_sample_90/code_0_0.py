# Coordinates of the points
D = (0, 3)
B = (0, 0)
C = (8, 0)

# Extracting coordinates
x1, y1 = D
x2, y2 = B
x3, y3 = C

# Calculating the area using the formula
area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

# Output the area
print(area)