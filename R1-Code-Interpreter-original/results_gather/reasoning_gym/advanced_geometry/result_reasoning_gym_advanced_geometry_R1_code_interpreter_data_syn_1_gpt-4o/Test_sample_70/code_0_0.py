from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A
eq1 = Eq(y - 1, (7/6)*(x - 9))

# Equation of altitude from B
eq2 = Eq(y, 7)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
print(solution)