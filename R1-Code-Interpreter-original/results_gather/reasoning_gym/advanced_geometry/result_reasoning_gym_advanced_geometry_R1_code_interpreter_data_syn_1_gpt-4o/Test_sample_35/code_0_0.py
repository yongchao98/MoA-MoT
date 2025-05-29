from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from C to AB (vertical line)
eq1 = Eq(x, -2)

# Equation of the altitude from A to BC
eq2 = Eq(y, 1 * (x + 7) + 3)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Print the solution
print((round(solution[x], 3), round(solution[y], 3)))