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

# Equations of the altitudes using point-slope form
# Altitude from A: y - 9 = slope_altitude_A * (x + 3)
# Altitude from B: y - 8 = slope_altitude_B * (x + 1)

# Solving the system of equations to find the intersection
# y = slope_altitude_A * (x + 3) + 9
# y = slope_altitude_B * (x + 1) + 8

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x + 3) + 9)
eq2 = Eq(y, slope_altitude_B * (x + 1) + 8)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
print(orthocenter)