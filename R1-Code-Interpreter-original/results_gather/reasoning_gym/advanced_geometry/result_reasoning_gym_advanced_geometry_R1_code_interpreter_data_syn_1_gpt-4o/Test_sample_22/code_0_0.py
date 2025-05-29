# Coordinates of the vertices
A = (4, 7)
B = (-1, 5)
C = (5, 2)

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
# y - y1 = m(x - x1) => y = mx + c
# Altitude from A
c_A = A[1] - slope_altitude_A * A[0]
# Altitude from B
c_B = B[1] - slope_altitude_B * B[0]

# Solve the system of equations to find the intersection (orthocenter)
# y = slope_altitude_A * x + c_A
# y = slope_altitude_B * x + c_B

# Set the equations equal to each other to find x
x_orthocenter = (c_B - c_A) / (slope_altitude_A - slope_altitude_B)
# Substitute x back into one of the equations to find y
y_orthocenter = slope_altitude_A * x_orthocenter + c_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))