from sympy import symbols, Eq, solve

# Define the variable
h = symbols('h')

# Define the equation based on similar triangles
equation = Eq(h / 3, (h - 1.5) / 0.5)

# Solve the equation for h
solution = solve(equation, h)

# Print the solution
print(solution[0])