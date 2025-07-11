# This script requires the 'spherogram' and 'sympy' libraries.
# You can install them using pip: pip install spherogram sympy

import spherogram
import sympy

# Define the braid beta in B_5
# sigma_i is represented by i, and sigma_i^-1 by -i.
# The number of strands is 5.
braid_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
b = spherogram.Braid(5, braid_word)

# Get the closure of the braid
K_beta_bar = b.closing()

# It's good practice to check that the closure is a knot (1 component)
if len(K_beta_bar.link_components) != 1:
    print("Warning: The braid closure is a link, not a knot.")
    # The Conway polynomial of a multi-component link is 0.
    poly_beta_bar_coeffs = [0]
else:
    # Compute the Conway polynomial for the braid closure
    poly_beta_bar_coeffs = K_beta_bar.conway_polynomial()

# Extract the z^2 coefficient for beta_bar
# The list is [c0, c1, c2, ...], so we need the element at index 2.
coeff_beta_bar = poly_beta_bar_coeffs[2] if len(poly_beta_bar_coeffs) > 2 else 0

# Get the knot 10_4 from the Rolfsen tables
K_10_4 = spherogram.Link('10_4')

# Compute the Conway polynomial for 10_4
poly_10_4_coeffs = K_10_4.conway_polynomial()

# Extract the z^2 coefficient for 10_4
coeff_10_4 = poly_10_4_coeffs[2] if len(poly_10_4_coeffs) > 2 else 0

# Calculate the difference
difference = coeff_beta_bar - coeff_10_4

# To display the full polynomials for clarity
z = sympy.Symbol('z')
poly_beta_bar_expr = sum(c * z**i for i, c in enumerate(poly_beta_bar_coeffs))
poly_10_4_expr = sum(c * z**i for i, c in enumerate(poly_10_4_coeffs))

print(f"The Alexander-Conway polynomial for the braid closure, Navla(beta_bar), is: {poly_beta_bar_expr}")
print(f"The Alexander-Conway polynomial for 10_4, Navla(10_4), is: {poly_10_4_expr}\n")

print("The problem is to find the difference in the z^2 coefficients.")
print(f"The z^2 coefficient of Navla(beta_bar) is: {coeff_beta_bar}")
print(f"The z^2 coefficient of Navla(10_4) is: {coeff_10_4}")
print(f"The difference is: {coeff_beta_bar} - {coeff_10_4} = {difference}")

<<< -6 >>>