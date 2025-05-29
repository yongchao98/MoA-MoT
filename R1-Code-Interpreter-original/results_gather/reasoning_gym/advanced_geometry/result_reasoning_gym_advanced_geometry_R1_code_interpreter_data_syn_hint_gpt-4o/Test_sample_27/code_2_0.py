from sympy import symbols, Eq, solve

# Define symbolic variables for the orthocenter coordinates
x, y = symbols('x y')

# Coordinates of the vertices
A = (-8, 10)
B = (4, -10)
C = (-6, -1)

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
# Altitude from A
eq1 = Eq(y, slope_altitude_A * (x - A[0]) + A[1])

# Altitude from B
eq2 = Eq(y, slope_altitude_B * (x - B[0]) + B[1])

# Solve the system of equations for the intersection of two altitudes
solution = solve((eq1, eq2), (x, y))

# Output the orthocenter coordinates
orthocenter_x = round(solution[x], 3)
orthocenter_y = round(solution[y], 3)
print(orthocenter_x, orthocenter_y)