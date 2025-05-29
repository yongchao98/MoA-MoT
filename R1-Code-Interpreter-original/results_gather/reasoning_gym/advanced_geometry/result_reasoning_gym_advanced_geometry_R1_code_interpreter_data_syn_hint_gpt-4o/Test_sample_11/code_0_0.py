from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A
eq1 = Eq(y + 6, (11/3)*(x + 1))

# Equation of the altitude from B
eq2 = Eq(y - 1, (3/5)*(x - 4))

# Solve the system of equations
orthocenter = solve((eq1, eq2), (x, y))
print(orthocenter)