from sympy import symbols, solve

# Define the symbol for the unknown radius
r3 = symbols('r3', positive=True)

# Define the curvatures
k1 = 1
k2 = 1/4
k3 = 1/r3
k4 = 0

# Descartes' Circle Theorem equation
equation = (k1 + k2 + k3 + k4)**2 - 2*(k1**2 + k2**2 + k3**2 + k4**2)

# Solve for r3
solution = solve(equation, r3)
print(solution)