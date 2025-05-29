from sympy import symbols, Eq, solve

# Define the variable
h = symbols('h')

# Equation from similar triangles
eq = Eq(h / 3, 1.5 / 2.5)

# Solve the equation
solution = solve(eq, h)

# Extract the height of the pole
pole_height = solution[0]
print(pole_height)