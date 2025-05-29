# Coordinates of the vertices
A = (9, 8)
B = (3, 5)
C = (-1, -6)

# Function to calculate slope
def slope(p1, p2):
    if p2[0] == p1[0]:  # vertical line
        return None
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = None if slope_BC is None else -1 / slope_BC
slope_altitude_B = None if slope_AC is None else -1 / slope_AC
slope_altitude_C = None if slope_AB is None else -1 / slope_AB

# Equations of the altitudes using point-slope form
# Altitude from A
if slope_altitude_A is not None:
    intercept_A = A[1] - slope_altitude_A * A[0]
else:
    intercept_A = A[0]  # vertical line

# Altitude from B
if slope_altitude_B is not None:
    intercept_B = B[1] - slope_altitude_B * B[0]
else:
    intercept_B = B[0]  # vertical line

# Solve the system of equations to find the intersection of two altitudes
# We will use the altitudes from A and B

# If both altitudes are not vertical
if slope_altitude_A is not None and slope_altitude_B is not None:
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
# If one of the altitudes is vertical
elif slope_altitude_A is None:
    x_orthocenter = intercept_A
    y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
elif slope_altitude_B is None:
    x_orthocenter = intercept_B
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))