from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equations of the lines
eq1 = Eq(y, -x + 4)
eq2 = Eq(y, -2*x + 12)

# Solve the system of equations
intersection = solve((eq1, eq2), (x, y))
intersection