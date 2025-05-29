# Coordinates of the vertices
A = (0, 5)
B = (-4, 2)
C = (3, 9)

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

# Equations of the altitudes
# Altitude from A (perpendicular to BC)
# y = -x + 5

# Altitude from B (perpendicular to AC)
# y = -3/4 * x + 2 + 4 * 3/4

# Solve the system of equations for the intersection of two altitudes
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, -x + 5)
eq2 = Eq(y, -3/4 * x + 2 + 4 * 3/4)

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)