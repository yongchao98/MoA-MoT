# Coordinates of the vertices
A = (1, 3)
B = (-7, 0)
C = (7, 6)

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

# Equations of the altitudes using point-slope form
# Altitude from A: y - 3 = slope_altitude_A * (x - 1)
# Altitude from B: y - 0 = slope_altitude_B * (x + 7)

# Solving the system of equations to find the intersection (orthocenter)
# y = slope_altitude_A * (x - 1) + 3
# y = slope_altitude_B * (x + 7)

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x - 1) + 3)
eq2 = Eq(y, slope_altitude_B * (x + 7))

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)