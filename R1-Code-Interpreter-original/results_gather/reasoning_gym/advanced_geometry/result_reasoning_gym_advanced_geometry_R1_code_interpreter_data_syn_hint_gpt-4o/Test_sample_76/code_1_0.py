from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A to BC
eq1 = Eq(y, -1/5 * x + 2/5)

# Equation of altitude from B to CA
eq2 = Eq(y, 2/5 * x + 14/5)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Round the solution to 3 decimal places
rounded_solution = {k: round(v, 3) for k, v in solution.items()}

# Print the rounded solution
print(rounded_solution)