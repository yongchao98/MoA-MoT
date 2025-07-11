import sympy as sp

# Define the symbols used in the equation
y, zeta_1, beta, k, H = sp.symbols('y zeta_1 beta k H')
psi_y = sp.Function('psi')(y)

# To represent the final equation, we construct it from its constituent parts,
# as derived from the physics of the problem.

# Part 1: The slip-dependent zeta potential at the bottom wall (y=-H/2)
zeta_slip_term = zeta_1 * (1 + beta * k)

# Part 2: The numerator of the hyperbolic ratio, based on the y-coordinate
numerator_term = sp.sinh(k * (H / 2 - y))

# Part 3: The denominator of the hyperbolic ratio, which normalizes the potential
denominator_term = sp.sinh(k * H)

# Assemble the final expression for the EDL potential psi(y)
potential_expression = zeta_slip_term * (numerator_term / denominator_term)

# Create the complete final equation psi(y) = expression
final_equation = sp.Eq(psi_y, potential_expression)

# Print out each part of the final equation and then the equation itself
print("The final expression for the Electrical Double-Layer (EDL) potential distribution psi(y) is constructed as follows:")
print("-" * 80)
print(f"Slip-dependent zeta potential term on the bottom wall: {zeta_slip_term}")
print(f"Spatial variation term (numerator): {numerator_term}")
print(f"Normalization term (denominator): {denominator_term}")
print("-" * 80)
print("The complete final equation is:")
# The sp.pprint function provides a nicely formatted output for the equation
sp.pprint(final_equation)