# The pyknotid library is required to run this code.
# You can install it by running: pip install pyknotid
import sympy
from pyknotid.representations import Braid
from pyknotid.spacecurves import Knot

# Step 1: Define the braid from the problem description
beta_generators = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
num_strands = 5
braid = Braid(beta_generators, num_strands=num_strands)

# Close the braid to get a knot
knot_beta_bar = Knot(braid.space_curve())

# Step 2: Calculate the Conway polynomial for the braid closure and get the z^2 coefficient
z = sympy.var('z')
poly_beta_bar = knot_beta_bar.conway_polynomial()
coeff_beta_bar_z2 = poly_beta_bar.coeff(z, 2)

# Step 3: Define the knot 10_4 and calculate its Conway polynomial
knot_10_4 = Knot(10, 4)
poly_10_4 = knot_10_4.conway_polynomial()
coeff_10_4_z2 = poly_10_4.coeff(z, 2)

# Step 4: Calculate and print the difference
difference = coeff_beta_bar_z2 - coeff_10_4_z2

print(f"The Alexander-Conway polynomial for the closure of beta is: {poly_beta_bar}")
print(f"The coefficient of z^2 for the closure of beta is: {coeff_beta_bar_z2}")
print("-" * 20)
print(f"The Alexander-Conway polynomial for 10_4 is: {poly_10_4}")
print(f"The coefficient of z^2 for 10_4 is: {coeff_10_4_z2}")
print("-" * 20)
print(f"The difference in the z^2 coefficients is: {coeff_beta_bar_z2} - ({coeff_10_4_z2}) = {difference}")

<<<-3>>>