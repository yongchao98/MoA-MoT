from sympy import symbols, Eq, solve

# Define the variables
a, b, c = symbols('a b c')

# Equations based on the problem
eq1 = Eq(a**2 + b**2, 3600)
eq2 = Eq((c-a)**2 + b**2, 6400)

# Solve the system of equations
solution = solve((eq1, eq2), (a, b, c))

# Extract the solution
a_val, b_val, c_val = solution[0]

# Calculate WZ
WZ = abs(c_val - a_val)

# Output the result
print(WZ)