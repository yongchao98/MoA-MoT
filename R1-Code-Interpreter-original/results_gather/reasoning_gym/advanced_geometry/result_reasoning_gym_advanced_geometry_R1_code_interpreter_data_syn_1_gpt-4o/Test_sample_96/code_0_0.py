# Coordinates of the vertices
A = (2, 6)
B = (1, -6)
C = (5, -10)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_AC
slope_altitude_C = -1 / slope_AB

# Equation of the altitude from A: y - y1 = m(x - x1)
# y = mx - mx1 + y1
# y = slope_altitude_A * x - slope_altitude_A * A[0] + A[1]
# Similarly for other altitudes

# Altitude from A
intercept_A = A[1] - slope_altitude_A * A[0]

# Altitude from B
intercept_B = B[1] - slope_altitude_B * B[0]

# Solve the system of equations for altitudes from A and B
# slope_altitude_A * x + intercept_A = slope_altitude_B * x + intercept_B
# (slope_altitude_A - slope_altitude_B) * x = intercept_B - intercept_A

x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))