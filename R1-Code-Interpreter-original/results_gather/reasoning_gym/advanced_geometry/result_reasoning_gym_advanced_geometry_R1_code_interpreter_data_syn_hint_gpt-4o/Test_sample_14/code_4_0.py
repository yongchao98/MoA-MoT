from sympy import symbols, Eq, solve

# Define symbolic variables
x, y = symbols('x y')

# Coordinates of the vertices
A = (2, -9)
B = (8, 8)
C = (-3, 10)

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

# Equations of the altitudes using point-slope form
# Altitude from A: y - y1 = m(x - x1)
eq_altitude_A = Eq(y, slope_altitude_A * (x - A[0]) + A[1])

# Altitude from B: y - y1 = m(x - x1)
eq_altitude_B = Eq(y, slope_altitude_B * (x - B[0]) + B[1])

# Solve for intersection of altitudes from A and B
solution = solve((eq_altitude_A, eq_altitude_B), (x, y))

# Extract and round the coordinates to 3 decimal places
x_orthocenter = round(solution[x], 3)
y_orthocenter = round(solution[y], 3)

print((x_orthocenter, y_orthocenter))