# Coordinates of the vertices
A = (-3, 9)
B = (-1, 8)
C = (-5, 1)

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
# y - 9 = slope_altitude_A * (x + 3)
# Equation of the altitude from B (y - y1 = m(x - x1))
# y - 8 = slope_altitude_B * (x + 1)

# Solving the system of equations to find the intersection
# y = slope_altitude_A * (x + 3) + 9
# y = slope_altitude_B * (x + 1) + 8

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x + 3) + 9)
eq2 = Eq(y, slope_altitude_B * (x + 1) + 8)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)