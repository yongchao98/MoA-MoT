import sympy as sp

# Define the symbols used in the problem's equations.
# K0: current sheet amplitude
# mu0: permeability of free space (air gap)
# mu: permeability of the magnetic material
# a: constant for spatial variation
# d: thickness of the air gap
# y: spatial coordinate
K0, mu0, mu, a, d, y = sp.symbols('K_0 mu_0 mu a d y')

# The final result is derived from magnetostatics principles. After solving
# Laplace's equation for the magnetic scalar potential and applying the
# boundary conditions, we find the magnetic field at the conductor's surface.
# The force per unit area is then found from the magnetic pressure P = B^2/(2*mu0),
# directed into the conductor (negative x-direction).

# We will construct the final expression piece by piece as requested.

# The numerator of the force magnitude expression
numerator = mu0 * K0**2 * sp.sin(a*y)**2

# The denominator of the force magnitude expression
# The base of the denominator is (cosh(ad) + (mu0/mu)*sinh(ad))
denominator_base = sp.cosh(a*d) + (mu0 / mu) * sp.sinh(a*d)
# The full denominator is 2 * (denominator_base)^2
denominator = 2 * denominator_base**2

# Now, we print the components of the final equation to show its structure.
print("The final equation for the force per unit area vector is composed of:")
print(f"Direction: Negative x-direction (-i_x)")
print(f"Numerator of magnitude: {numerator}")
print(f"Denominator of magnitude: {denominator}")

print("\nCombining these parts, the full symbolic expression for the force per unit area is:")

# We use sp.pprint for a clean, formatted mathematical representation.
# First, display the magnitude.
force_magnitude = numerator / denominator
print("\nMagnitude:")
sp.pprint(force_magnitude)

# Then, display the full vector expression.
print("\nFull vector (f/area):")
# Create a symbolic representation of the unit vector for printing
i_x = sp.Symbol('i_x')
sp.pprint(-force_magnitude * i_x, use_unicode=True)