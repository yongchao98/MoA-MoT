from sympy import symbols, Eq, solve

# Define symbolic variables for the orthocenter coordinates
x, y = symbols('x y')

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

# Equations of the altitudes
# Altitude from A
if slope_altitude_A is not None:
    eq1 = Eq(y - A[1], slope_altitude_A * (x - A[0]))
else:
    eq1 = Eq(x, A[0])

# Altitude from B
if slope_altitude_B is not None:
    eq2 = Eq(y - B[1], slope_altitude_B * (x - B[0]))
else:
    eq2 = Eq(x, B[0])

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Extract and round the solution
x_orthocenter = round(solution[x], 3)
y_orthocenter = round(solution[y], 3)

print((x_orthocenter, y_orthocenter))