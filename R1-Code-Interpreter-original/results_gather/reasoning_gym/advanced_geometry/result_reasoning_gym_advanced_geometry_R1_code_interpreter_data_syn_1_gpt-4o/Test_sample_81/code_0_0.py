# Coordinates of the vertices
A = (10, -2)
B = (-10, 5)
C = (7, 1)

# Function to calculate the slope of a line given two points
def slope(p1, p2):
    if p2[0] - p1[0] == 0:
        return None  # Vertical line
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC if slope_BC is not None else None
slope_altitude_B = -1 / slope_AC if slope_AC is not None else None
slope_altitude_C = -1 / slope_AB if slope_AB is not None else None

# Function to find the y-intercept of a line given a point and slope
def y_intercept(point, slope):
    if slope is None:
        return None
    return point[1] - slope * point[0]

# Equations of the altitudes
intercept_A = y_intercept(A, slope_altitude_A)
intercept_B = y_intercept(B, slope_altitude_B)

# Solve the system of equations for the intersection of two altitudes
# y = slope_altitude_A * x + intercept_A
# y = slope_altitude_B * x + intercept_B

# If both slopes are not None, solve for x
if slope_altitude_A is not None and slope_altitude_B is not None:
    x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
    y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A
else:
    # Handle vertical line cases
    if slope_altitude_A is None:
        x_orthocenter = A[0]
        y_orthocenter = slope_altitude_B * x_orthocenter + intercept_B
    elif slope_altitude_B is None:
        x_orthocenter = B[0]
        y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))