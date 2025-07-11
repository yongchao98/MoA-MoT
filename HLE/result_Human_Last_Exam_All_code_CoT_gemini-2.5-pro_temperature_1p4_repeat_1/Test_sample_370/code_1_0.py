import math

# The final result for the total cross section (sigma) is derived from the principles
# of quantum field theory, assuming a high-energy collision where particle masses
# can be neglected. The result depends on the interaction coupling constant 'g'
# and the energy of a single incoming fermion 'E' in the center-of-mass frame.

# The derived formula is: sigma = (3 * g^4) / (128 * pi * E^2)

# Below we define the numerical constants present in the final equation.
numerator_constant = 3
denominator_constant = 128

# The symbolic parts are 'g' (coupling constant), 'pi', and 'E' (energy).

# The code will now print each part of the equation as requested.
print("The final equation for the total cross section sigma is constructed as follows:")
print("-" * 60)
print(f"Numerator of the equation:")
print(f"  - A numerical constant: {numerator_constant}")
print(f"  - A symbolic part: g**4 (the coupling constant to the fourth power)")

print(f"\nDenominator of the equation:")
print(f"  - A numerical constant: {denominator_constant}")
print(f"  - Symbolic parts: pi * E**2 (pi times the energy squared)")

print("-" * 60)
print("\nPutting it all together, the final expression is:")
print(f"sigma = ({numerator_constant} * g**4) / ({denominator_constant} * pi * E**2)")
