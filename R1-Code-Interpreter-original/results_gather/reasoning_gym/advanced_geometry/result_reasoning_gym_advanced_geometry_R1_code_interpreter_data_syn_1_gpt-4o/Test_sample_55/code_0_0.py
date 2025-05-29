# Coordinates of the vertices
A = (-4, 7)
B = (-7, -2)
C = (0, 10)

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
# Altitude from A: y - y1 = m(x - x1)
# Altitude from B: y - y1 = m(x - x1)

# Altitude from A
# y - 7 = slope_altitude_A * (x + 4)
# y = slope_altitude_A * x + slope_altitude_A * 4 + 7

# Altitude from B
# y + 2 = slope_altitude_B * (x + 7)
# y = slope_altitude_B * x + slope_altitude_B * 7 - 2

# Solving the system of equations
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * x + slope_altitude_A * 4 + 7)
eq2 = Eq(y, slope_altitude_B * x + slope_altitude_B * 7 - 2)

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)