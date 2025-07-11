import math

def solve_for_m():
    """
    This function calculates the probability Pm for a given integer m.
    A sample value m=10 is used here. You can change this value.
    """
    # You can change the value of m here
    m = 10

    # Total number of pairs (i, j) that can be chosen from 4m+2 elements
    # which is C(4m+2, 2) = (4m+2)(4m+1)/2 = (2m+1)(4m+1)
    total_pairs = (2 * m + 1) * (4 * m + 1)

    # Calculate the number of "good" pairs based on the parity of m.
    
    if m % 2 != 0:
        # If m is odd, only partitions where sub-sequences have a common difference of 1 are considered.
        # The number of such pairs is C(m+2, 2).
        good_pairs = math.comb(m + 2, 2)
    else:
        # If m is even, we consider partitions with delta=1 and delta=2.
        # Let m = 2k
        k = m // 2
        
        # Count for pairs allowing delta=1 partitions:
        c1_count = math.comb(m + 2, 2)
        
        # Count for pairs allowing delta=2 partitions:
        c2_count = (k + 1)**2
        
        # Count the intersection of the two sets of pairs to avoid double counting:
        c1_c2_intersect = math.comb(k + 2, 2)
        
        # Total good pairs is the size of the union of the two sets.
        good_pairs = c1_count + c2_count - c1_c2_intersect

    # As requested, printing the components of the final probability equation for the given m.
    print(f"For m = {m}, the probability P_m is given by the fraction:")
    print(f"Numerator (number of valid pairs): {good_pairs}")
    print(f"Denominator (total number of pairs): {total_pairs}")
    print(f"So, P_{m} = {good_pairs} / {total_pairs}")

solve_for_m()