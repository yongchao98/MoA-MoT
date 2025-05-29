from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A
eq1 = Eq(y - 7, -1/5 * (x - 6))

# Equation of the altitude from B (vertical line)
eq2 = Eq(x, -5)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
print(solution)