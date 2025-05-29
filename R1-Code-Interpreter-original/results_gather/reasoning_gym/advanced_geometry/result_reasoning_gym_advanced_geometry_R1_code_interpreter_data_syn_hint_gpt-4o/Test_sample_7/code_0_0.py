from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A
eq1 = Eq(y - 2, -7/6 * (x + 2))

# Equation of altitude from B (vertical line)
eq2 = Eq(x, -10)

# Solve for intersection of eq1 and eq2
intersection1 = solve((eq1, eq2), (x, y))

# Equation of altitude from C
eq3 = Eq(y - 2, -4/3 * (x + 3))

# Solve for intersection of eq1 and eq3
intersection2 = solve((eq1, eq3), (x, y))

# Print the results
print(intersection1)
print(intersection2)