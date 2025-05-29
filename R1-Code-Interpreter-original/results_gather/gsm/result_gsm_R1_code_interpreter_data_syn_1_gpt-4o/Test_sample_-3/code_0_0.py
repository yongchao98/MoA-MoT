from sympy import symbols, Eq, solve

# Define the variables
c, w = symbols('c w')

# Define the equations
eq1 = Eq(c + w, 3016226)
eq2 = Eq(2*c + 4*w, 70)

# Solve the system of equations
solution = solve((eq1, eq2), (c, w))
print(solution)