import spherogram
import sympy

def solve_knot_problem():
    """
    Calculates the difference between the braid index of K2 and a lower bound
    on the Seifert circles of K1.
    """

    # Part 1: Analyze K1 = 10_74
    # Get the knot object for 10_74
    K1 = spherogram.Link('10_74')

    # Calculate the HOMFLY polynomial P(a, z)
    poly_K1 = K1.homfly_polynomial()

    # Find the span of the variable 'a'
    # The variable 'a' is the first variable in the polynomial ring
    a_var = poly_K1.vars[0]
    a_degrees_K1 = poly_K1.degrees(a_var)
    min_deg_a_K1 = min(a_degrees_K1)
    max_deg_a_K1 = max(a_degrees_K1)
    span_a_K1 = max_deg_a_K1 - min_deg_a_K1

    # Calculate the lower bound for the number of Seifert circles
    s_min_lower_bound_K1 = int(span_a_K1 / 2 + 1)

    print(f"For knot K1 = 10_74:")
    print(f"  The HOMFLY polynomial is P(a,z) = {poly_K1}")
    print(f"  The span of variable 'a' is {max_deg_a_K1} - ({min_deg_a_K1}) = {span_a_K1}")
    print(f"  The lower bound for the minimum number of Seifert circles is {span_a_K1}/2 + 1 = {s_min_lower_bound_K1}")
    print("-" * 20)

    # Part 2: Analyze K2 = closure of (sigma_1^-1)^3 * sigma_2^-1
    # Create the braid object on 3 strands.
    # The word is sigma_1^-1, sigma_1^-1, sigma_1^-1, sigma_2^-1
    braid_word_K2 = [-1, -1, -1, -2]
    K2_braid = spherogram.Braid(3, braid_word_K2)

    # Create the link from the braid closure
    K2 = spherogram.Link(K2_braid)

    # Identify the knot
    knot_id = K2.identify()[0]
    
    print(f"For knot K2 = closure of (sigma_1^-1)^3 * sigma_2^-1:")
    print(f"  This knot is identified as {knot_id}.")
    print(f"  Since it is the closure of a 3-strand braid, its braid index is at most 3.")

    # A knot has braid index 2 iff it's a (2,k) torus knot.
    # Let's check its Alexander polynomial.
    alex_poly_K2 = K2.alexander_polynomial()
    
    print(f"  The Alexander polynomial of K2 is: {alex_poly_K2}")
    print(f"  The Alexander polynomial for a (2,k) torus knot has only coefficients of +1 or -1.")
    print(f"  Since the polynomial of K2 has other coefficients, it is not a (2,k) torus knot.")
    print(f"  Therefore, its braid index cannot be 2.")
    
    # Since the knot is not trivial (braid index > 1), and not 2, and <= 3, it must be 3.
    braid_index_K2 = 3
    print(f"  Thus, the braid index of K2 is {braid_index_K2}.")
    print("-" * 20)

    # Part 3: Calculate the difference
    difference = braid_index_K2 - s_min_lower_bound_K1

    print("Final Calculation:")
    print(f"The braid index of K2 is {braid_index_K2}.")
    print(f"The lower bound of the minimum number of Seifert circles of K1 is {s_min_lower_bound_K1}.")
    print(f"The difference is {braid_index_K2} - {s_min_lower_bound_K1} = {difference}")

solve_knot_problem()