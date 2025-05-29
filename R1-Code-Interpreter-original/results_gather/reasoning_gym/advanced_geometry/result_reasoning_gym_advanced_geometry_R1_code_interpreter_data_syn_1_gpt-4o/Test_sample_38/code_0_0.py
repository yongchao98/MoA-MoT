# Coordinates of the vertices
A = (-1, 4)
B = (4, -6)
C = (0, -7)

# Function to calculate the slope
def slope(p1, p2):
    if p2[0] - p1[0] == 0:
        return None  # vertical line
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC if slope_BC is not None else 0
slope_altitude_B = -1 / slope_AC if slope_AC is not None else 0
slope_altitude_C = -1 / slope_AB if slope_AB is not None else 0

# Equations of the altitudes
# Altitude from A: y - y1 = m(x - x1)
# Altitude from B: y - y2 = m(x - x2)
# Altitude from C: y - y3 = m(x - x3)

# Altitude from A
x1, y1 = A
if slope_altitude_A is not None:
    intercept_A = y1 - slope_altitude_A * x1
else:
    intercept_A = x1

# Altitude from B
x2, y2 = B
if slope_altitude_B is not None:
    intercept_B = y2 - slope_altitude_B * x2
else:
    intercept_B = x2

# Altitude from C
x3, y3 = C
if slope_altitude_C is not None:
    intercept_C = y3 - slope_altitude_C * x3
else:
    intercept_C = x3

# Solve the system of equations to find the orthocenter
# We can solve for the intersection of two altitudes, say from A and B

if slope_altitude_A is not None and slope_altitude_B is not None:
    # y = slope_altitude_A * x + intercept_A
    # y = slope_altitude_B * x + intercept_B
    # Set the equations equal to each other to find x
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
elif slope_altitude_A is None:
    # Altitude from A is vertical
    x_orthocenter = intercept_A
    y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
elif slope_altitude_B is None:
    # Altitude from B is vertical
    x_orthocenter = intercept_B
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))