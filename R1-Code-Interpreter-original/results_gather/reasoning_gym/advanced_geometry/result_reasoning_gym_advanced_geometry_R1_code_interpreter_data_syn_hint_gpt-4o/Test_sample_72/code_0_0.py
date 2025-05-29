from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equations of the altitudes
eq1 = Eq(y, -5)  # Altitude from A
eq2 = Eq(y, x - 1)  # Altitude from B

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
print(solution)