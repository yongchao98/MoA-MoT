from sympy import symbols, Eq, solve

# Define the variable
n = symbols('n', integer=True, positive=True)

# Define the equation
equation = Eq((1 + n) * (1/2)**n, 3/16)

# Solve the equation
solution = solve(equation, n)

# Print the solution
print(solution)