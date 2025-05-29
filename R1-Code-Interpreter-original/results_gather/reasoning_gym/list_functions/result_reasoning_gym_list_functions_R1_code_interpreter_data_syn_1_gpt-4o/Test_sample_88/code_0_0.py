from sympy import symbols, Eq, solve

# Define the variables
a, b = symbols('a b')

# Define the equations based on the examples
eq1 = Eq(a * 25 + b, 46)
eq2 = Eq(a * 29 + b, 54)

# Solve the equations
solution = solve((eq1, eq2), (a, b))
print(solution)