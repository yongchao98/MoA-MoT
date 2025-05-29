from sympy import symbols, Eq, solve

# Define the variable
n = symbols('n', integer=True, positive=True)

# Define the equation based on the given condition
equation = Eq(n - (n * (n - 3)) / 2, 0)

# Solve the equation
solution = solve(equation, n)

# Print the solution
print(solution)