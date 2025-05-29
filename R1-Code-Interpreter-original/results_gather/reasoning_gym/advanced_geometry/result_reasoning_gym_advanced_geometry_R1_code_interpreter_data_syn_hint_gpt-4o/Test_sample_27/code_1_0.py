# Coordinates of the vertices
A = (-8, 10)
B = (4, -10)
C = (-6, -1)

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

# Equations of the altitudes
# Altitude from A: y - y1 = m(x - x1)
# Altitude from B: y - y2 = m(x - x2)
# Altitude from C: y - y3 = m(x - x3)

# Altitude from A
intercept_A = A[1] - slope_altitude_A * A[0]

# Altitude from B
intercept_B = B[1] - slope_altitude_B * B[0]

# Solve the system of equations for the intersection of two altitudes
# y = slope_altitude_A * x + intercept_A
# y = slope_altitude_B * x + intercept_B

# Solving for x
x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Output the orthocenter coordinates
print(round(x_orthocenter, 3), round(y_orthocenter, 3))