import spherogram
import sympy

def solve_knot_problem():
    """
    This script solves the knot theory problem by performing the following steps:
    1. Identifies the knot K2 from its braid representation and determines its braid index.
    2. Calculates the lower bound for the minimum number of Seifert circles of knot K1 (10_74)
       from its HOMFLY polynomial.
    3. Computes the difference between these two values.
    """
    print("--- Part 1: Analyzing Knot K2 ---")
    
    # K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1 in B_3.
    # The braid word is [-1, -1, -1, -2] on 3 strands.
    braid_K2 = spherogram.Braid(3, [-1, -1, -1, -2])
    K2_link = spherogram.Link(braid_K2)

    # Identify the knot using spherogram.
    # The identify() method returns a list; we take the first and best identification.
    K2 = K2_link.identify()[0]
    print(f"The knot K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1, which is identified as the {K2.name()} knot.")

    # Determine the braid index of K2.
    # The given representation is on 3 strands, so the braid index b(K2) <= 3.
    # Knots with braid index 2 are torus knots, which have cyclotomic Alexander polynomials.
    # The Alexander polynomial of 5_2 is 2*t^2 - 3*t + 2, which is not cyclotomic.
    # Thus, b(K2) > 2. Since b(K2) <= 3, the braid index must be 3.
    braid_index_K2 = 3
    print(f"The braid index of K2 ({K2.name()}) is {braid_index_K2}.\n")

    print("--- Part 2: Analyzing Knot K1 ---")

    # K1 is the 10_74 knot.
    K1 = spherogram.Knot('10_74')
    print(f"The knot K1 is the {K1.name()} knot.")

    # Calculate the HOMFLY polynomial for K1.
    homfly_poly_K1 = K1.homfly_polynomial()
    print(f"The HOMFLY polynomial P(v, z) of K1 is: {homfly_poly_K1}")

    # The lower bound for the number of Seifert circles s(K) is given by
    # s(K) >= v-span(P(K)) + 1. We need to find the v-span.
    v = homfly_poly_K1.gens[0]
    p_poly = sympy.poly(homfly_poly_K1, v)
    
    # Extract degrees of v.
    coeffs = p_poly.as_dict()
    v_degrees = [d[0] for d in coeffs.keys()]
    max_deg_v = max(v_degrees)
    min_deg_v = min(v_degrees)
    v_span = max_deg_v - min_deg_v
    
    print(f"The maximum power of v in the polynomial is {max_deg_v}.")
    print(f"The minimum power of v in the polynomial is {min_deg_v}.")
    
    # Calculate the lower bound for the number of Seifert circles.
    seifert_circles_lower_bound_K1 = v_span + 1
    print(f"The v-span of the polynomial is {max_deg_v} - {min_deg_v} = {v_span}.")
    print(f"The lower bound of the minimum number of Seifert circles of K1 is {v_span} + 1 = {seifert_circles_lower_bound_K1}.\n")

    print("--- Part 3: Calculating the Final Difference ---")

    # Calculate the difference between the two quantities.
    difference = abs(braid_index_K2 - seifert_circles_lower_bound_K1)
    
    print(f"The value for K2 (braid index) is: {braid_index_K2}")
    print(f"The value for K1 (Seifert circles lower bound) is: {seifert_circles_lower_bound_K1}")
    print("The final equation for the difference is:")
    print(f"{braid_index_K2} - {seifert_circles_lower_bound_K1} = {braid_index_K2 - seifert_circles_lower_bound_K1}")
    print(f"The absolute difference is {difference}.")
    
    return difference

if __name__ == '__main__':
    solve_knot_problem()
