# Coordinates of the vertices
A = (9, 8)
B = (3, 5)
C = (-1, -6)

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
x1, y1 = A
m1 = slope_altitude_A
# y - 8 = m1(x - 9)
# y = m1*x - 9*m1 + 8

# Altitude from B
x2, y2 = B
m2 = slope_altitude_B
# y - 5 = m2(x - 3)
# y = m2*x - 3*m2 + 5

# Solve the system of equations to find the intersection
# m1*x - 9*m1 + 8 = m2*x - 3*m2 + 5

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(m1*x - 9*m1 + 8, y)
eq2 = Eq(m2*x - 3*m2 + 5, y)

# Solve for x and y
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)