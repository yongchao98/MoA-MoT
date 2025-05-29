# Coordinates of the vertices
A = (-1, -6)
B = (4, 1)
C = (-7, 4)

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

# Function to calculate the y-intercept of a line given a point and slope
def y_intercept(point, slope):
    return point[1] - slope * point[0]

# Equations of the altitudes
# y = slope * x + intercept
intercept_A = y_intercept(A, slope_altitude_A)
intercept_B = y_intercept(B, slope_altitude_B)
intercept_C = y_intercept(C, slope_altitude_C)

# Solve the system of equations to find the intersection of two altitudes
# Altitude from A and B
# y = slope_altitude_A * x + intercept_A
# y = slope_altitude_B * x + intercept_B

# Set the equations equal to each other and solve for x
x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the results to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))