from sympy import symbols, Eq, solve

# Define the variables
a, d = symbols('a d')

# Define the equations
eq1 = Eq(a + 3*d, 120)
eq2 = Eq(4*a + 6*d, 360)

# Solve the system of equations
solution = solve((eq1, eq2), (a, d))

# Extract the smallest angle
smallest_angle = solution[a]
print(smallest_angle)