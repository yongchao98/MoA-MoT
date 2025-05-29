from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equations of the altitudes
eq1 = Eq(y + 6, (6/7)*(x + 4))
eq2 = Eq(y + 3, (9/4)*(x + 1))

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))

# Extract the coordinates and round them
x_orthocenter = round(solution[x], 3)
y_orthocenter = round(solution[y], 3)

print((x_orthocenter, y_orthocenter))