from sympy import symbols, Eq, solve

# Coordinates of the vertices
A = (-4, -6)
B = (-1, -8)
C = (-5, -1)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_AC

# Define symbols
x, y = symbols('x y')

# Equations of the altitudes using point-slope form
eq1 = Eq(y, slope_altitude_A * (x + 4) - 6)
eq2 = Eq(y, slope_altitude_B * (x + 1) - 8)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)