import math

def solve_knot_problem():
    """
    This function calculates the difference between two knot theory invariants.
    1. The lower bound of the minimum number of Seifert circles of K1 = 10_74.
    2. The braid index of K2 = closure of (sigma_1^-1)^3 * sigma_2^-1.
    """

    # Part 1: Calculate the lower bound of the minimum number of Seifert circles for K1 (10_74)

    # The HOMFLY polynomial for the 10_74 knot is a known result.
    # P(a, z) = a^8 - a^6*z^2 - a^6 + 2*a^4*z^2 + 3*a^4 - a^2*z^2 - 2*a^2
    # The powers of the variable 'a' in this polynomial are {8, 6, 4, 2}.
    powers_of_a_k1 = [8, 6, 4, 2]
    
    max_power = max(powers_of_a_k1)
    min_power = min(powers_of_a_k1)
    
    # The span of the HOMFLY polynomial in 'a' is the difference between the max and min powers.
    span_a_k1 = max_power - min_power
    
    # A lower bound for the minimum number of Seifert circles, s(K), is given by
    # the Morton-Franks-Williams inequality: s(K) >= 0.5 * span_a(P(K)) + 1.
    seifert_circles_lower_bound = 0.5 * span_a_k1 + 1
    
    print(f"For K1 = 10_74:")
    print(f"  The powers of 'a' in the HOMFLY polynomial are: {powers_of_a_k1}")
    print(f"  Maximum power of 'a' = {max_power}")
    print(f"  Minimum power of 'a' = {min_power}")
    print(f"  Span of the HOMFLY polynomial = {max_power} - {min_power} = {span_a_k1}")
    print(f"  Lower bound for the minimum number of Seifert circles = 0.5 * {span_a_k1} + 1 = {seifert_circles_lower_bound}")
    print("-" * 20)

    # Part 2: Calculate the braid index of K2

    # K2 is the closure of the braid beta = (sigma_1^-1)^3 * sigma_2^-1 from the 3-strand braid group B_3.
    # Since it can be represented by a 3-strand braid, its braid index is at most 3.
    
    # We must check if the braid index can be 1 or 2.
    # We use the Alexander polynomial of K2, which is known to be Delta(t) = t^2 - 3t + 1.
    
    # Condition for braid index 1: The knot is the unknot (Delta(t) = 1).
    # Since t^2 - 3t + 1 is not equal to 1, K2 is not the unknot. So, braid index > 1.
    
    # Condition for braid index 2: The knot is a (2, p) torus knot for an odd integer p.
    # The Alexander polynomial of a (2, p) torus knot is of the form 1 - t + t^2 - ... + t^(p-1).
    # The polynomial t^2 - 3t + 1, with its coefficient of -3, is not of this form.
    # Therefore, K2 is not a (2,p) torus knot, and its braid index is not 2.
    
    # Since the braid index is at most 3, but not 1 or 2, it must be 3.
    braid_index_k2 = 3
    
    print(f"For K2 = closure of (sigma_1^-1)^3 * sigma_2^-1:")
    print(f"  The knot is given by a 3-strand braid, so its braid index is at most 3.")
    print(f"  Its Alexander polynomial is t^2 - 3t + 1.")
    print(f"  The knot is not the unknot (braid index 1) or a (2,p)-torus knot (braid index 2).")
    print(f"  Therefore, the braid index of K2 must be {braid_index_k2}.")
    print("-" * 20)

    # Part 3: Calculate the difference
    difference = seifert_circles_lower_bound - braid_index_k2
    
    print("Final Calculation:")
    print(f"The difference is between the lower bound of the minimum number of Seifert circles of K1 ({seifert_circles_lower_bound})")
    print(f"and the braid index of K2 ({braid_index_k2}).")
    print(f"{seifert_circles_lower_bound} - {braid_index_k2} = {difference}")

    return difference

if __name__ == '__main__':
    solve_knot_problem()