# Coordinates of the vertices
A = (8, -7)
B = (0, -3)
C = (10, 9)

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

# Equation of the altitude from A (y - y1 = m(x - x1))
# y + 7 = slope_altitude_A * (x - 8)
# y = slope_altitude_A * x - slope_altitude_A * 8 - 7

# Equation of the altitude from B (y - y1 = m(x - x1))
# y + 3 = slope_altitude_B * (x - 0)
# y = slope_altitude_B * x - 3

# Solving the system of equations for altitudes from A and B
# slope_altitude_A * x - y = slope_altitude_A * 8 + 7
# slope_altitude_B * x - y = 3

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(slope_altitude_A * x - y, slope_altitude_A * 8 + 7)
eq2 = Eq(slope_altitude_B * x - y, 3)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)