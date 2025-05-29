from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A
eq1 = Eq(y + 5, -9/2 * (x - 8))

# Equation of the altitude from B
eq2 = Eq(y - 2, 4/9 * (x + 5))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)