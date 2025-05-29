# Coordinates of the vertices
A = (10, -2)
B = (-10, 5)
C = (7, 1)

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
# Altitude from A: y - A[1] = slope_altitude_A * (x - A[0])
# Altitude from B: y - B[1] = slope_altitude_B * (x - B[0])

# Solving the system of equations to find the intersection (orthocenter)
# y = slope_altitude_A * (x - A[0]) + A[1]
# y = slope_altitude_B * (x - B[0]) + B[1]

# Let's solve these equations
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, slope_altitude_A * (x - A[0]) + A[1])
eq2 = Eq(y, slope_altitude_B * (x - B[0]) + B[1])

# Solve for x and y
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)