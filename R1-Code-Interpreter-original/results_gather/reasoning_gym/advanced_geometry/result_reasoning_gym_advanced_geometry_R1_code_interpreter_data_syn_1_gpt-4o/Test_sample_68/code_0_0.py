# Coordinates of the vertices
A = (4, 4)
B = (-5, 2)
C = (1, 3)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slope of BC
slope_BC = slope(B, C)

# Slope of altitude from A (perpendicular to BC)
slope_altitude_A = -1 / slope_BC

# Equation of altitude from A: y - y1 = m(x - x1)
# y - 4 = slope_altitude_A * (x - 4)
# y = slope_altitude_A * x - slope_altitude_A * 4 + 4

# Slope of AC
slope_AC = slope(A, C)

# Slope of altitude from B (perpendicular to AC)
slope_altitude_B = -1 / slope_AC

# Equation of altitude from B: y - y1 = m(x - x1)
# y - 2 = slope_altitude_B * (x + 5)
# y = slope_altitude_B * x + slope_altitude_B * 5 + 2

# Solve the system of equations to find the intersection
# y = slope_altitude_A * x - slope_altitude_A * 4 + 4
# y = slope_altitude_B * x + slope_altitude_B * 5 + 2

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * x - slope_altitude_A * 4 + 4)
eq2 = Eq(y, slope_altitude_B * x + slope_altitude_B * 5 + 2)

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)