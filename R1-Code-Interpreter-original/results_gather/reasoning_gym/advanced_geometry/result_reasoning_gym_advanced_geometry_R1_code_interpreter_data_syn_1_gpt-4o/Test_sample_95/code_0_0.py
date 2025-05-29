from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A to BC
eq1 = Eq(y - 8, -6 * (x - 9))

# Equation of altitude from B to AC
eq2 = Eq(y + 7, -19/18 * (x - 8))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)