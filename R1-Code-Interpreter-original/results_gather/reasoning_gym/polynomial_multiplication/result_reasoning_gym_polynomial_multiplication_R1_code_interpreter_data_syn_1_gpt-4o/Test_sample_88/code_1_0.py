from sympy import symbols, expand

# Define the symbol
z = symbols('z')

# Define the polynomials
first_poly = 28*z**2 - 82*z
second_poly = -9*z**3 - 25

# Perform the multiplication
simplified_expression = expand(first_poly * second_poly)

# Print the simplified expression
print(simplified_expression)