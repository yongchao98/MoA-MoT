# Define the variable
r = symbols('r')

# Define the revised equation
revised_equation = Eq(3*r + (r + 2), 42)

# Solve the revised equation
red_chickens_revised = solve(revised_equation, r)

# Print the result
print(red_chickens_revised)