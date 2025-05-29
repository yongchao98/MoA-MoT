# Coordinates of the vertices
A = (3, 1)
B = (-9, 8)
C = (-10, 7)

# Function to calculate the slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
m_BC = slope(B, C)
m_CA = slope(C, A)
m_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
m_altitude_A = -1 / m_BC
m_altitude_B = -1 / m_CA
m_altitude_C = -1 / m_AB

# Equations of the altitudes using point-slope form
# Altitude from A: y - 1 = m_altitude_A * (x - 3)
# Altitude from B: y - 8 = m_altitude_B * (x + 9)

# Solving the system of equations to find the intersection
# y = m_altitude_A * (x - 3) + 1
# y = m_altitude_B * (x + 9) + 8

from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, m_altitude_A * (x - 3) + 1)
eq2 = Eq(y, m_altitude_B * (x + 9) + 8)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)