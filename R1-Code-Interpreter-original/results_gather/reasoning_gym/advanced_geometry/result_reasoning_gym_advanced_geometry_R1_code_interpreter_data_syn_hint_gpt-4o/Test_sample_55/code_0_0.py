from sympy import symbols, Eq, solve

# Define variables
x, y = symbols('x y')

# Equation of the altitude from A: y - 7 = (-7/12)(x + 4)
eq1 = Eq(y - 7, (-7/12)*(x + 4))

# Equation of the altitude from B: y + 2 = (-4/3)(x + 7)
eq2 = Eq(y + 2, (-4/3)*(x + 7))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Print the solution
print(solution)