# Coordinates of the vertices
A = (2, 6)
B = (1, -6)
C = (5, -10)

# Function to calculate the slope of a line given two points
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

# Function to find the y-intercept of a line given a point and a slope
def y_intercept(point, slope):
    return point[1] - slope * point[0]

# Equations of the altitudes
# y = slope * x + intercept
intercept_A = y_intercept(A, slope_altitude_A)
intercept_B = y_intercept(B, slope_altitude_B)

# Solve the system of equations for the altitudes from A and B
# slope_altitude_A * x + intercept_A = slope_altitude_B * x + intercept_B
# (slope_altitude_A - slope_altitude_B) * x = intercept_B - intercept_A
x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))