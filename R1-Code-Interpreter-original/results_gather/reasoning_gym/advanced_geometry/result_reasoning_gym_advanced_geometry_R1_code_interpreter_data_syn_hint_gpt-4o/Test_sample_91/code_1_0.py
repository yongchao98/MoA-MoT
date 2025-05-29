from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A
eq1 = Eq(y + 8, -10/7 * (x + 2))

# Equation of the altitude from B
eq2 = Eq(y + 2, -5 * (x - 3))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Round the solution to three decimal places
orthocenter = {key: round(value, 3) for key, value in solution.items()}

# Print the rounded solution
print(orthocenter)