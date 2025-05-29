# Define the points
A = (9, -2)
B = (-6, 1)
C = (-8, -6)

# Function to calculate slope
def slope(p1, p2):
    if p1[0] == p2[0]:  # vertical line
        return None
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (perpendicular slopes)
slope_altitude_A = -1 / slope_BC if slope_BC is not None else 0
slope_altitude_B = -1 / slope_AC if slope_AC is not None else 0
slope_altitude_C = -1 / slope_AB if slope_AB is not None else 0

# Equations of the altitudes using point-slope form: y - y1 = m(x - x1)
# Altitude from A
x1, y1 = A
if slope_altitude_A is not None:
    intercept_A = y1 - slope_altitude_A * x1
else:
    intercept_A = x1  # vertical line

# Altitude from B
x2, y2 = B
if slope_altitude_B is not None:
    intercept_B = y2 - slope_altitude_B * x2
else:
    intercept_B = x2  # vertical line

# Solve for intersection of two altitudes
# y = slope_altitude_A * x + intercept_A
# y = slope_altitude_B * x + intercept_B

if slope_altitude_A is not None and slope_altitude_B is not None:
    # Solve for x
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    # Solve for y
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
elif slope_altitude_A is None:  # Altitude from A is vertical
    x_orthocenter = intercept_A
    y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
elif slope_altitude_B is None:  # Altitude from B is vertical
    x_orthocenter = intercept_B
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))