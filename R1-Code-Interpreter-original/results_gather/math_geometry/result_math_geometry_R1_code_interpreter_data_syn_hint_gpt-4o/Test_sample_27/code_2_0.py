from sympy import symbols, sqrt, solve

# Define the variables
R, h, a, c = symbols('R h a c')

# Given values
R_value = 12
h_value = 24

# Expression for the radius of the inscribed sphere
r_expr = (R * h) / (sqrt(R**2 + h**2) + R)

# Substitute the given values
r_value = r_expr.subs({R: R_value, h: h_value})

# Set up the equation r = a(sqrt(c) - 1)
eq = a * (sqrt(c) - 1) - r_value

# Solve for a and c
solution = solve(eq, (a, c))

# Calculate a + c
result = solution[0][0] + solution[0][1]

print(result)