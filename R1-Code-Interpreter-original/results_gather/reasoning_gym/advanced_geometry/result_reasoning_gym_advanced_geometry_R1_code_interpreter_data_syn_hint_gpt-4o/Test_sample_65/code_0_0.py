# Coordinates of the vertices
A = (3, 1)
B = (-9, 8)
C = (-10, 7)

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
# Altitude from A: y - 1 = slope_altitude_A * (x - 3)
# Altitude from B: y - 8 = slope_altitude_B * (x + 9)

# Solving the system of equations to find the intersection
# y = slope_altitude_A * (x - 3) + 1
# y = slope_altitude_B * (x + 9) + 8

# Equating the two expressions for y
# slope_altitude_A * (x - 3) + 1 = slope_altitude_B * (x + 9) + 8

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x - 3) + 1)
eq2 = Eq(y, slope_altitude_B * (x + 9) + 8)

# Solve the equations
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)