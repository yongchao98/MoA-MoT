from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A: y = -x + 5
eq1 = Eq(y, -x + 5)

# Corrected equation of the altitude from B: y = -3/4 * x - 1
eq2 = Eq(y, -3/4 * x - 1)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)