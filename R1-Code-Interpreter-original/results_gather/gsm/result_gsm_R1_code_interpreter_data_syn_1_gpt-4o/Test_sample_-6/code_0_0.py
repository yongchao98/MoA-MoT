from sympy import symbols, Eq, solve

# Define the variable
r = symbols('r')

# Define the equation
equation = Eq(3*r + 1801110*(r + 2), 42)

# Solve the equation
solution = solve(equation, r)

# Print the solution
print(solution)