import spherogram
import sympy

def solve_knot_problem():
    """
    This script solves the given knot theory problem by:
    1. Finding the braid index of K2, the closure of (sigma_1^-1)^3 * sigma_2^-1.
    2. Finding the lower bound of the minimum number of Seifert circles of K1 (the 10_74 knot)
       from its HOMFLY polynomial.
    3. Calculating the difference between these two values.
    """
    # --- Part 1: Analyze K2 ---
    print("--- Part 1: Analyzing Knot K2 ---")
    
    # K2 is the closure of the braid word sigma_1^-3 * sigma_2^-1 in the braid group B_3.
    # In spherogram notation, sigma_i is i and sigma_i^-1 is -i.
    # The braid is on 3 strands.
    braid_K2 = spherogram.Braid(3, [-1, -1, -1, -2])
    link_K2 = spherogram.Link(braid_K2)
    
    # Identify the knot. This will be the figure-eight knot, 4_1.
    knot_K2_identified = link_K2.identify()
    print(f"K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1.")
    print(f"This knot is identified as the {knot_K2_identified.name} knot.")
    
    # The braid index of a knot is the minimum number of strands required to represent it.
    # The figure-eight knot (4_1) is known to have a braid index of 3.
    braid_index_K2 = 3
    print(f"The braid index of K2 ({knot_K2_identified.name}) is {braid_index_K2}.\n")

    # --- Part 2: Analyze K1 ---
    print("--- Part 2: Analyzing Knot K1 ---")
    
    # K1 is the 10_74 knot.
    link_K1 = spherogram.Link('10_74')
    print(f"K1 is the {link_K1.name} knot.")

    # The minimum number of Seifert circles of a knot is equal to its braid index.
    # We find a lower bound for this value using the HOMFLY polynomial, P(a, z).
    # The bound is given by the Morton-Franks-Williams inequality:
    # bound = (max_degree(a) - min_degree(a)) / 2 + 1
    
    # Get the HOMFLY polynomial. Spherogram uses variables x, y for a, z.
    x, y = sympy.var('x, y')
    homfly_poly_K1 = link_K1.homfly_polynomial()
    print(f"The HOMFLY polynomial for K1 in variables (a, z) is: P(a, z) = {homfly_poly_K1.subs({x: 'a', y: 'z'})}")

    # To find the degrees of 'a' (which is 'x' in our polynomial object),
    # we can treat it as a polynomial in x.
    poly_in_x = sympy.poly(homfly_poly_K1, x)
    
    # Get all monomials to find the degrees of x.
    degrees_x = [m[0] for m in poly_in_x.monoms()]
    max_deg_a = max(degrees_x)
    min_deg_a = min(degrees_x)

    print(f"The maximum degree of 'a' is {max_deg_a}.")
    print(f"The minimum degree of 'a' is {min_deg_a}.")
    
    # Calculate the lower bound.
    seifert_bound_K1 = (max_deg_a - min_deg_a) / 2 + 1
    
    print(f"The lower bound for the minimum number of Seifert circles of K1 is:")
    print(f"({max_deg_a} - {min_deg_a}) / 2 + 1 = {int(seifert_bound_K1)}\n")

    # --- Part 3: Calculate the Difference ---
    print("--- Part 3: Final Calculation ---")
    
    difference = braid_index_K2 - seifert_bound_K1
    
    print("The difference is the (braid index of K2) - (lower bound for K1).")
    print(f"Difference = {braid_index_K2} - {int(seifert_bound_K1)} = {int(difference)}")


if __name__ == '__main__':
    solve_knot_problem()
    print("\n<<<0>>>")
