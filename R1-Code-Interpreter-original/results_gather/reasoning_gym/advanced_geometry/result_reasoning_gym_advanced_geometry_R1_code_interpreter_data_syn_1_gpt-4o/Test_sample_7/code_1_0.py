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
slope_altitude_A = None if slope_BC is None else (-1 / slope_BC if slope_BC != 0 else None)
slope_altitude_B = None if slope_CA is None else (-1 / slope_CA if slope_CA != 0 else None)
slope_altitude_C = None if slope_AB is None else (-1 / slope_AB if slope_AB != 0 else None)

# Equation of altitude from A (using point-slope form)
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
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
elif slope_altitude_A is None:
    x_orthocenter = intercept_A
    y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
else:
    x_orthocenter = intercept_B
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))