import sympy
from pyknotid.catalogue import get_knot
from pyknotid.catalogue.identify import from_braid

# Step 1: Define the braid word for beta
# sigma_i is represented by i, and sigma_i^-1 by -i.
beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

# Step 2: Create knot objects
# Create a knot from the closure of the braid beta
k_beta = from_braid(beta_word)

# Create a knot object for 10_4 from the catalogue
k_10_4 = get_knot('10_4')

# Step 3: Calculate the Alexander-Conway polynomials
poly_beta = k_beta.alexander_conway_polynomial()
poly_10_4 = k_10_4.alexander_conway_polynomial()

# Step 4: Extract the z^2 coefficients
# Get the symbolic variable 'z' from the polynomial
try:
    z = list(poly_beta.free_symbols)[0]
except IndexError: # This handles constant polynomials
    z = sympy.Symbol('z')

coeff_beta_z2 = poly_beta.coeff(z**2)
coeff_10_4_z2 = poly_10_4.coeff(z**2)

# Step 5: Calculate the difference and print the result
difference = coeff_beta_z2 - coeff_10_4_z2

print(f"The Alexander-Conway polynomial for the braid closure (beta_bar) is: {poly_beta}")
print(f"The Alexander-Conway polynomial for the knot 10_4 is: {poly_10_4}")
print("-" * 20)
print(f"The z^2 coefficient for beta_bar is: {coeff_beta_z2}")
print(f"The z^2 coefficient for 10_4 is: {coeff_10_4_z2}")
print("-" * 20)
print("The difference between the coefficients is:")
# The final equation with the numbers, as requested
print(f"{coeff_beta_z2} - ({coeff_10_4_z2}) = {difference}")
