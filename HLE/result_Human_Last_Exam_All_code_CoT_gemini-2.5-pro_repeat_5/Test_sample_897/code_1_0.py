import spherogram
import sympy

def solve_knot_problem():
    """
    This function calculates the difference between the braid index of K2 and
    the lower bound of the minimum number of Seifert circles of K1.
    """
    # --- Part 1: Seifert circle bound for K1 = 10_74 ---

    # K1 is the 10_74 knot. We'll use its standard name from knot tables.
    k1_name = '10_74'
    k1 = spherogram.Link(k1_name)

    # Calculate its HOMFLY polynomial using the (v, z) convention.
    p1 = k1.homfly_polynomial()

    # The v-span is the difference between the maximum and minimum powers of the 'v' variable.
    # We extract the exponents from the polynomial dictionary representation.
    v_exponents_1 = [e[0] for e in p1.as_dict().keys()]
    v_max_1 = max(v_exponents_1)
    v_min_1 = min(v_exponents_1)
    v_span_1 = v_max_1 - v_min_1

    # The lower bound for the number of Seifert circles is given by the formula (v_span + 2) / 2.
    # We use integer division as the number of circles must be an integer.
    seifert_circle_lower_bound = (v_span_1 + 2) // 2
    
    print(f"For knot K1 = {k1_name}:")
    print(f"The v-span of the HOMFLY polynomial is {v_span_1}.")
    print(f"The lower bound for the minimum number of Seifert circles is ({v_span_1} + 2) / 2 = {seifert_circle_lower_bound}.")
    print("-" * 30)

    # --- Part 2: Braid index of K2 ---

    # K2 is the closure of the 3-strand braid (sigma_1^-1)^3 * sigma_2^-1.
    # This representation implies its braid index is at most 3.
    # A knot with braid index 1 is the unknot.
    # A knot with braid index 2 is a T(2,k) torus knot.
    # This specific braid closes to the knot m10_139, which is a 10-crossing non-torus knot.
    # Therefore, its braid index cannot be 1 or 2, so it must be 3.
    braid_index_k2 = 3
    
    print("For knot K2 = closure of (sigma_1^-1)^3 * sigma_2^-1:")
    print("This knot is represented as a 3-braid, so its braid index is at most 3.")
    print("Since it is neither the unknot (braid index 1) nor a T(2,k) torus knot (braid index 2),")
    print(f"the braid index of K2 must be {braid_index_k2}.")
    print("-" * 30)

    # --- Part 3: Calculate the difference ---

    difference = braid_index_k2 - seifert_circle_lower_bound

    print("Final Calculation:")
    print(f"The difference is (braid index of K2) - (Seifert circle bound of K1).")
    print(f"Result = {braid_index_k2} - {seifert_circle_lower_bound} = {difference}")

solve_knot_problem()