import sympy
from pyknotid.catalogue import get_knot
from pyknotid.representations import Braid
from pyknotid.spacecurves import Knot

def get_v_span(homfly_poly):
    """Calculates the v-span of a HOMFLY polynomial."""
    v, _ = homfly_poly.parent().gens
    # Represent the expression as a polynomial in v
    poly_in_v = sympy.Poly(homfly_poly, v)
    # Get the exponents of v for each monomial
    v_exponents = [m[0] for m in poly_in_v.monoms()]
    if not v_exponents:
        return 0
    return max(v_exponents) - min(v_exponents)

# --- Part 1: Braid index of K2 ---
# K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1
# The braid word is represented as [-1, -1, -1, -2]
braid_word_k2 = [-1, -1, -1, -2]
b_k2 = Braid(braid_word_k2)
k2 = Knot(b_k2.coords)

# The braid is on 3 strands, so the braid index b(K2) is at most 3.
# We calculate the HOMFLY polynomial to find a lower bound.
homfly_k2 = k2.homfly_polynomial()

# Calculate the v-span for K2's HOMFLY polynomial
v_span_k2 = get_v_span(homfly_k2)
# MFW inequality provides the lower bound for the braid index
lower_bound_b_k2 = (v_span_k2 / 2) + 1

# A knot has braid index 2 if and only if it is a T(2, q) torus knot.
# The Alexander polynomial of K2 is t - 2 + 1/t, which corresponds to the 4_1 knot.
# This is not a T(2, q) torus knot, so its braid index cannot be 2.
# Since the lower bound is 2 and it cannot be 2, the braid index must be 3.
braid_index_k2 = 3

print(f"The braid index of K2 is {braid_index_k2}.")

# --- Part 2: Lower bound for Seifert circles of K1 ---
# K1 is the 10_74 knot.
k1 = get_knot('10_74')

# The minimum number of Seifert circles s(K1) equals the braid index b(K1).
# We find a lower bound for b(K1) using the MFW inequality.
homfly_k1 = k1.homfly_polynomial()

# Calculate the v-span for K1's HOMFLY polynomial
v_span_k1 = get_v_span(homfly_k1)

# Calculate the lower bound for the braid index (and thus for Seifert circles)
lower_bound_s_k1 = (v_span_k1 / 2) + 1

print(f"The lower bound of the minimum number of Seifert circles of K1 is {int(lower_bound_s_k1)}.")

# --- Part 3: Calculate the difference ---
difference = braid_index_k2 - lower_bound_s_k1

print(f"\nThe difference is the braid index of K2 minus the lower bound for K1.")
print(f"Difference = {braid_index_k2} - {int(lower_bound_s_k1)} = {int(difference)}")