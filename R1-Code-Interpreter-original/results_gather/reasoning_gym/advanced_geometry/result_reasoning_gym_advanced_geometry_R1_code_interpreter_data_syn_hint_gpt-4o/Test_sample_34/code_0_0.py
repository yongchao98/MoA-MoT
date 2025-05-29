from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equations of the altitudes
eq1 = Eq(y, 2*x)  # Altitude from C to AB
eq2 = Eq(y, (1/3)*x + 10/3)  # Altitude from A to BC
eq3 = Eq(y, 4)  # Altitude from B to CA

# Solve the system of equations
intersection1 = solve((eq1, eq2), (x, y))
intersection2 = solve((eq1, eq3), (x, y))
intersection3 = solve((eq2, eq3), (x, y))

# Print the intersection points
print(intersection1)
print(intersection2)
print(intersection3)