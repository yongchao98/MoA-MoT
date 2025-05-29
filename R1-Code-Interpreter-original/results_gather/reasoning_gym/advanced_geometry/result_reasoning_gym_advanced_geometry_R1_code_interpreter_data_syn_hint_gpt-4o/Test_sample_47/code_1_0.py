# Re-evaluating the orthocenter calculation

# Slopes of the altitudes
slope_altitude_A = 6
slope_altitude_B = -7/4

# Equations of the altitudes
# Altitude from A: y - 1 = 6(x - 4)
# Altitude from B: y + 5 = -7/4(x - 9)

# Solving the system of equations
from sympy import symbols, Eq, solve

x, y = symbols('x y')
eq1 = Eq(y, 6 * (x - 4) + 1)
eq2 = Eq(y, -7/4 * (x - 9) - 5)

# Solve the equations
orthocenter = solve((eq1, eq2), (x, y))
orthocenter_rounded = (round(orthocenter[x], 3), round(orthocenter[y], 3))

print(orthocenter_rounded)