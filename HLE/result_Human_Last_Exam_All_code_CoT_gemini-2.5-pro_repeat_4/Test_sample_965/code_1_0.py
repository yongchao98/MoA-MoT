import sympy

# Define the symbols
g, h, gamma_c, pi = sympy.symbols('g h gamma_c pi')

# The rate for making a photon is given by the formula:
# Rate = 8 * pi * g^2 / (h * gamma_c)
# We will print out the components of this formula.

numerator_coeff = 8
numerator_vars = [pi, g**2]

denominator_vars = [h, gamma_c]

# Using print to display the final equation structure
print("The rate for making a photon is given by the formula:")
print("Rate = (numerator) / (denominator)")
print(f"numerator = {numerator_coeff} * {pi} * {g**2}")
print(f"denominator = {h} * {gamma_c}")

# This corresponds to the expression 8 * pi * g^2 / (h * gamma_c)
# Let's print the final assembled expression string as well
final_expression_str = f"{numerator_coeff} * pi * g^2 / (h * gamma_c)"
print("\nFinal Expression:")
print(final_expression_str)