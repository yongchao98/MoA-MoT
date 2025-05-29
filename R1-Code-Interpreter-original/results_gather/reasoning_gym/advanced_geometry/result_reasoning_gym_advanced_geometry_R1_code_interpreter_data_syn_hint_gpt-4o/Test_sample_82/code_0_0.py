# Coordinates of the vertices
A = (8, -7)
B = (0, -3)
C = (10, 9)

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
# Altitude from A: y - y1 = m(x - x1)
# Altitude from B: y - y1 = m(x - x1)
# We will find the intersection of altitudes from A and B

# Altitude from A: y + 7 = slope_altitude_A * (x - 8)
# Altitude from B: y + 3 = slope_altitude_B * (x - 0)

# Solving the system of equations
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y + 7, slope_altitude_A * (x - 8))
eq2 = Eq(y + 3, slope_altitude_B * (x - 0))

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)