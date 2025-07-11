import sympy

# Define symbols for the chemical potentials and charge
mu_2 = sympy.Symbol('μ_2')
mu_3 = sympy.Symbol('μ_3')
e = sympy.Symbol('e')

# The formula for the second voltage plateau is modeled as the average of the
# potentials associated with Stage 2 and Stage 3.
# V = - (mu_avg) / e
# mu_avg = (μ_2 + μ_3) / 2
# So, V = -(μ_2 + μ_3) / (2*e)

# Construct the expression
voltage_expression = -(mu_2 + mu_3) / (2 * e)

# Print the formula in a human-readable format.
# The problem asks to output each number in the final equation.
# The numbers are the indices '2' and '3', and the coefficient '2'.
print(f"The formula is: {voltage_expression}")
print("\nOr, in a more standard text format:")
print("-(μ_2 + μ_3) / (2 * e)")
