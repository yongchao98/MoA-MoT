from sympy import symbols, solve

# Define the variable
r3 = symbols('r3', positive=True)

# Curvatures
k1 = 1
k2 = 0.25
k3 = 1/r3

# Descartes' Circle Theorem equation
equation = (k1 + k2 + k3)**2 - 2*(k1**2 + k2**2 + k3**2)

# Solve for r3
solution = solve(equation, r3)
print(solution)