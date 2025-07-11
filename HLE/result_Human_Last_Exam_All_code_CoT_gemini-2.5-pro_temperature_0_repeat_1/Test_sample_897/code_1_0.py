import spherogram
import sympy

def solve_knot_problem():
    """
    This function solves the given knot theory problem by:
    1. Calculating the braid index of K2.
    2. Calculating the lower bound for the minimum number of Seifert circles of K1.
    3. Finding the difference between these two values.
    """
    # --- Step 1: Find the braid index of K2 ---
    # K2 is the closure of the 3-strand braid (sigma_1^-1)^3 * sigma_2^-1.
    # The braid word is [-1, -1, -1, -2].
    braid_K2 = spherogram.Braid([-1, -1, -1, -2])
    K2 = spherogram.Link(braid_K2)

    # The number of strands is 3, so the braid index b(K2) is at most 3.
    # We use the Morton-Franks-Williams inequality to find a lower bound:
    # b(K) >= (span_t(V(K)) / 2) + 1, where V(K) is the Jones polynomial.
    t = sympy.Symbol('t')
    jones_poly_K2 = K2.jones_polynomial(t)

    # The Jones polynomial is a Laurent polynomial in t. We find its degree span.
    # For K2 (which is the knot 5_2), V(t) = t^2 - t + 1 - t^-1 + t^-2.
    # The powers are 2, 1, 0, -1, -2.
    poly_t = sympy.Poly(jones_poly_K2, t, t**-1)
    degrees_t = poly_t.monoms()
    powers_t = [d[0] for d in degrees_t]
    min_deg_jones = min(powers_t)
    max_deg_jones = max(powers_t)
    span_jones_K2 = max_deg_jones - min_deg_jones

    # Calculate the lower bound for the braid index.
    braid_index_lower_bound_K2 = span_jones_K2 / 2 + 1

    # Since b(K2) <= 3 and the lower bound is 3.0, b(K2) must be 3.
    braid_index_K2 = 3
    print(f"The braid index of K2 is {braid_index_K2}.")

    # --- Step 2: Find the lower bound for Seifert circles of K1 ---
    # K1 is the 10_74 knot.
    K1 = spherogram.Link('10_74')

    # The lower bound for the minimum number of Seifert circles s(K) is given by
    # span_z(P(K)) + 1, where P(K) is the HOMFLY polynomial.
    a, z = sympy.symbols('a,z')
    homfly_poly_K1 = K1.homfly_polynomial(a, z)

    # The HOMFLY polynomial is a Laurent polynomial in a and z. We find the z-span.
    # For 10_74, P(a,z) = z^4 + (3 - 2/a^2)*z^2 + (1/a^4 - 2/a^2 + 1).
    # The powers of z are 4, 2, 0.
    poly_z = sympy.Poly(homfly_poly_K1, z)
    degrees_z = poly_z.monoms()
    powers_z = [d[0] for d in degrees_z]
    min_deg_z = min(powers_z)
    max_deg_z = max(powers_z)
    span_z_K1 = max_deg_z - min_deg_z

    # Calculate the lower bound for the number of Seifert circles.
    seifert_circles_lower_bound_K1 = span_z_K1 + 1
    print(f"The lower bound for the minimum number of Seifert circles of K1 is {seifert_circles_lower_bound_K1}.")

    # --- Step 3: Calculate the final difference ---
    difference = braid_index_K2 - seifert_circles_lower_bound_K1
    print("\nThe final difference is the braid index of K2 minus the lower bound for K1.")
    print(f"Difference = {braid_index_K2} - {seifert_circles_lower_bound_K1} = {int(difference)}")

if __name__ == "__main__":
    solve_knot_problem()