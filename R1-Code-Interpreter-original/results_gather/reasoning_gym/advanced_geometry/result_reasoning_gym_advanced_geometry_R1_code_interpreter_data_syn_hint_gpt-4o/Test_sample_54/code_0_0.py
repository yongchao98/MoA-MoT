from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Define the equations of the altitudes
eq1 = Eq(y, -8/9 * x + 56/9)
eq2 = Eq(y, 1/10 * x - 7)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Print the solution rounded to 3 decimal places
orthocenter = (round(solution[x], 3), round(solution[y], 3))
print(orthocenter)