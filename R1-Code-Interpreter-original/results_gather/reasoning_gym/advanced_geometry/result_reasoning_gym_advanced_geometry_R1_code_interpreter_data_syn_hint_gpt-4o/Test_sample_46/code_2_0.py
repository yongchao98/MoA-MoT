from sympy import symbols, Eq, solve

# Define the symbols
x, y = symbols('x y')

# Coordinates of the vertices
A = (-7, -1)
B = (-2, 3)
C = (-9, 10)

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
eq1 = Eq(y - A[1], slope_altitude_A * (x - A[0]))
eq2 = Eq(y - B[1], slope_altitude_B * (x - B[0]))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Extract the coordinates of the orthocenter
x_orthocenter = round(solution[x], 3)
y_orthocenter = round(solution[y], 3)

print((x_orthocenter, y_orthocenter))