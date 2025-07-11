import sympy

# Define the symbols for the equation
g, E, pi = sympy.symbols('g E pi')
sigma = sympy.Symbol('sigma')

# Coefficients and powers from the derived formula
numerator_coeff = 3
g_power = 4
denominator_coeff = 128
E_power = 2

# Construct the formula for the total cross section
# σ = (3 * g^4) / (128 * π * E^2)
cross_section_formula = (numerator_coeff * g**g_power) / (denominator_coeff * pi * E**E_power)

# Create a formatted string for the equation
equation_str = f"{sigma} = {sympy.pretty(cross_section_formula, use_unicode=False)}"

# Print the final equation, showing each number as requested
print("The total cross section σ for the scattering of two fermions in the high-energy limit is:")
print(equation_str)