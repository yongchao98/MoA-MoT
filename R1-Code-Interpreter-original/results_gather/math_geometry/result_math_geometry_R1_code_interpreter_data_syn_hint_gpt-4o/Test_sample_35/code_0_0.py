# Vertices of the triangle
x1, y1 = 0, -4
x2, y2 = 0, 16
x3, y3 = 4, 4

# Calculate the area using the formula
area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
print(area)