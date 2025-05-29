from sympy import symbols, Eq, solve

# Define the variables
a, b = symbols('a b')

# Define the equations based on the examples
eq1 = Eq(a * 97 + b, 190)
eq2 = Eq(a * 44 + b, 84)
eq3 = Eq(a * 63 + b, 122)
eq4 = Eq(a * 24 + b, 44)

# Solve the system of equations
solution = solve((eq1, eq2), (a, b))

# Print the solution
print(solution)