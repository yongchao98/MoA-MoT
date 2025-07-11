import sympy

# Define the symbols
g, h, gamma_c, pi = sympy.symbols('g h gamma_c pi')

# The formula for the photon production rate based on a common result in cavity QED
# that matches one of the options.
# Γ = 8 * π * g^2 / (h * γ_c)
rate_numerator = 8 * pi * g**2
rate_denominator = h * gamma_c
rate = rate_numerator / rate_denominator

# To present it clearly, let's print the components
print("The rate (Γ) for making a photon is given by the formula:")
print("Γ = (8 * π * g^2) / (h * γ_c)")
print("\nWhich corresponds to the expression:")
# Use sympy.pretty_print for a nice mathematical layout
sympy.pretty_print(rate)
