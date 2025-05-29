# Coordinates of the vertices
A = (-1, 4)
B = (4, -6)
C = (0, -7)

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
# y - 4 = slope_altitude_A * (x + 1)
# Equation of the altitude from B (y + 6 = slope_altitude_B * (x - 4))
# Solve these two equations to find the orthocenter

# Solving the system of equations
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y - 4, slope_altitude_A * (x + 1))
eq2 = Eq(y + 6, slope_altitude_B * (x - 4))

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)