from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A to BC
eq1 = Eq(y - 2, -1/5 * (x + 8))

# Equation of the altitude from B to CA
eq2 = Eq(y - 2, 0.4 * (x + 2))

# Solve the system of equations
intersection = solve((eq1, eq2), (x, y))

# Print the intersection point
print(intersection)