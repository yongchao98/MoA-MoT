import numpy as np
from pyknotid.representations import Braid
from pyknotid.catalogue import get_knot

def get_z2_coeff(poly):
    """
    Extracts the coefficient of the z^2 term from a numpy.poly1d object
    representing a Conway polynomial.
    
    A polynomial p(z) = a_n*z^n + ... + a_0 is represented by poly.coeffs
    as [a_n, ..., a_0]. The coefficient of z^k is at index (poly.order - k).
    """
    order = poly.order
    # If the polynomial degree is less than 2, the z^2 coefficient is 0.
    if order < 2:
        return 0
    
    # Get the coefficient of z^2 (where k=2)
    coeff = poly.coeffs[order - 2]
    # Conway coefficients are integers.
    return int(coeff)

# Step 1: Define the braid beta from the problem description.
# In pyknotid, the generator sigma_i is represented by index i-1.
# The braid word is s4^-1 s4^-1 s3^-1 s4 s3^-1 s2 s1^-1 s3^-1 s2^-1 s2^-1 s2^-1 s1^-1
braid_word = [
    (3, -1), (3, -1), (2, -1), (3, 1), (2, -1), (1, 1), (0, -1),
    (2, -1), (1, -1), (1, -1), (1, -1), (0, -1)
]
beta = Braid(braid_word, num_strands=5)

# Step 2: Create the knot k_beta from the closure of the braid beta.
k_beta = beta.get_knot()

# Step 3: Compute the Conway polynomial for the closure of beta.
nabla_beta = k_beta.conway_polynomial()

# Step 4: Extract the z^2 coefficient for the knot from beta.
c2_beta = get_z2_coeff(nabla_beta)

# Step 5: Get the knot 10_4 from the pyknotid catalogue.
k_10_4 = get_knot('10_4')

# Step 6: Compute the Conway polynomial for 10_4.
nabla_10_4 = k_10_4.conway_polynomial()

# Step 7: Extract the z^2 coefficient for 10_4.
c2_10_4 = get_z2_coeff(nabla_10_4)

# Step 8: Calculate the difference and print the results with the final equation.
difference = c2_beta - c2_10_4

print(f"The Alexander-Conway polynomial for the closure of beta is: {nabla_beta}")
print(f"The z^2 coefficient for the closure of beta is: {c2_beta}")
print(f"\nThe Alexander-Conway polynomial for 10_4 is: {nabla_10_4}")
print(f"The z^2 coefficient for 10_4 is: {c2_10_4}")
print(f"\nThe difference in the z^2 coefficients is: {c2_beta} - {c2_10_4} = {difference}")
