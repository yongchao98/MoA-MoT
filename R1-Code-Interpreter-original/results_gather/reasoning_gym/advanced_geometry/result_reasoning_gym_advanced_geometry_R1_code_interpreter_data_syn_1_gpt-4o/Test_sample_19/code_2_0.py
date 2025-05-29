from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of altitude from C to AB
eq3 = Eq(y - 7, -11/10 * (x + 3))

# Equation of altitude from B to AC (vertical line)
eq2 = Eq(x, -5)

# Solve the system of equations
solution2 = solve((eq3, eq2), (x, y))
print(solution2)