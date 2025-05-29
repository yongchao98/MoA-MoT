# Coordinates of the vertices
A = (-2, 2)
B = (-10, -4)
C = (-3, 2)

# Function to calculate slope
def slope(p1, p2):
    if p2[0] == p1[0]:  # vertical line
        return None
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_AB = slope(A, B)
slope_BC = slope(B, C)
slope_CA = slope(C, A)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = None if slope_BC is None else -1 / slope_BC
slope_altitude_B = None if slope_CA is None else -1 / slope_CA
slope_altitude_C = None if slope_AB is None else -1 / slope_AB

# Equation of altitude from A (using point-slope form)
# y - y1 = m(x - x1) => y = m*x + (y1 - m*x1)
if slope_altitude_A is not None:
    intercept_A = A[1] - slope_altitude_A * A[0]
else:
    intercept_A = A[0]  # x = constant

# Equation of altitude from B
if slope_altitude_B is not None:
    intercept_B = B[1] - slope_altitude_B * B[0]
else:
    intercept_B = B[0]  # x = constant

# Solve for intersection of altitudes from A and B
if slope_altitude_A is not None and slope_altitude_B is not None:
    # y = m1*x + c1 and y = m2*x + c2
    # m1*x + c1 = m2*x + c2 => (m1 - m2)*x = c2 - c1
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
elif slope_altitude_A is None:
    # Altitude from A is vertical, x = intercept_A
    x_orthocenter = intercept_A
    y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
else:
    # Altitude from B is vertical, x = intercept_B
    x_orthocenter = intercept_B
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))