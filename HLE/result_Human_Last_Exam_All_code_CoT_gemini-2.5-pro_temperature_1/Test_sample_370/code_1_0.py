import sympy

# Define the symbols used in the final equation
g, E, pi = sympy.symbols('g E pi')

# From the derivation, the total cross section sigma is calculated.
# We will construct and print the final expression.

# The coefficient in the denominator is 32
coeff = 32

# Construct the numerator and denominator strings for clear printing
numerator_str = "g**4"
denominator_str = f"{coeff} * pi * E**2"

print("The total cross section for the scattering of two fermions in the lowest order is:")
print(f"σ = {numerator_str} / ({denominator_str})")
print("\nIn this equation:")
print(f"  - g is the coupling constant of the interaction.")
print(f"  - E is the energy of each fermion in the center-of-mass frame.")
print(f"  - The number in the denominator is {coeff}.")