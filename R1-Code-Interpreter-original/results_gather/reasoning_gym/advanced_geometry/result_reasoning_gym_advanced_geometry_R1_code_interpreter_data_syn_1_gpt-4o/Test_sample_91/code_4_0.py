# Coordinates of the vertices
A = (-2, -8)
B = (3, -2)
C = (-7, -9)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_AC

# Equations of the altitudes using point-slope form
# Altitude from A: y - y1 = m(x - x1)
# Altitude from B: y - y1 = m(x - x1)

# Altitude from A
x1, y1 = A
m1 = slope_altitude_A
# y - y1 = m1(x - x1) => y = m1*x - m1*x1 + y1
b1 = -m1 * x1 + y1

# Altitude from B
x2, y2 = B
m2 = slope_altitude_B
# y - y2 = m2(x - x2) => y = m2*x - m2*x2 + y2
b2 = -m2 * x2 + y2

# Solve for intersection (orthocenter)
# m1*x + b1 = m2*x + b2
# (m1 - m2)x = b2 - b1
x_orthocenter = (b2 - b1) / (m1 - m2)
y_orthocenter = m1 * x_orthocenter + b1

# Round to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))