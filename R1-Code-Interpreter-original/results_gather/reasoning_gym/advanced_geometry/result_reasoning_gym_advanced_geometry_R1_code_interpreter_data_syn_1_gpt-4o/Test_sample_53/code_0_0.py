from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from A
eq1 = Eq(y + 7, (3/8)*(x - 9))

# Equation of altitude from B
eq2 = Eq(y + 4, (12/11)*x)

# Solve the system of equations
solution = solve((eq1, eq2), (x, y))
orthocenter = (round(solution[x], 3), round(solution[y], 3))

print(orthocenter)