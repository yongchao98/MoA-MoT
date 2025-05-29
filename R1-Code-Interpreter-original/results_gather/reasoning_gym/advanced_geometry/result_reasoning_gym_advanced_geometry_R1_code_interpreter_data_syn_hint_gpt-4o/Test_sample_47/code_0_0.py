# Coordinates of the vertices
A = (4, 1)
B = (9, -5)
C = (-3, -3)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_AB = slope(A, B)
slope_BC = slope(B, C)
slope_CA = slope(C, A)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_CA
slope_altitude_C = -1 / slope_AB

# Equations of the altitudes using point-slope form
# Altitude from A: y - 1 = slope_altitude_A * (x - 4)
# Altitude from B: y + 5 = slope_altitude_B * (x - 9)

# Solving the system of equations to find the intersection (orthocenter)
# y = slope_altitude_A * (x - 4) + 1
# y = slope_altitude_B * (x - 9) - 5

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x - 4) + 1)
eq2 = Eq(y, slope_altitude_B * (x - 9) - 5)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)