from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A to BC
eq1 = Eq(y - 3, 1 * (x + 7))

# Equation of the altitude from B to CA
eq2 = Eq(y - 3, -5/6 * (x - 4))

# Solve for the intersection of the altitudes
intersection = solve((eq1, Eq(x, -2)), (x, y))
print(intersection)