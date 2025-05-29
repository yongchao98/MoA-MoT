from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from A
eq1 = Eq(y, -8/9 * x + 56/9)

# Equation of the altitude from B
eq2 = Eq(y, 1/10 * x - 9)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Round the solution to 3 decimal places
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)