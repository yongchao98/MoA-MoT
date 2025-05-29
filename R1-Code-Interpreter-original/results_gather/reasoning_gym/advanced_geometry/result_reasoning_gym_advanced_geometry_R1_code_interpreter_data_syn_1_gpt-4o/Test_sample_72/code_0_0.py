from sympy import symbols, Eq, solve

# Define variables
x, y = symbols('x y')

# Equations of the altitudes
eq1 = Eq(y, -5)  # Altitude from A to BC
eq2 = Eq(y + 6, x + 7)  # Altitude from B to AC
eq3 = Eq(y + 8, 3 * (x + 7))  # Altitude from C to AB

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
print(solution)