import math

def solve():
    """
    Calculates the number of sequences and the parameters a1, b1 for the largest term.
    """

    # Part 1: Count the number of sequences.
    # This is interpreted as the popcount of N = 10^100.
    # popcount(10^100) = popcount( (2*5)^100 ) = popcount( 2^100 * 5^100 )
    # Since multiplying by a power of 2 only shifts the bits, it doesn't change the popcount.
    # So, popcount(10^100) = popcount(5^100).
    
    num = 5**100
    
    # In Python 3.10+, we can use num.bit_count()
    # For compatibility, we can use bin(num).count('1')
    try:
        count = num.bit_count()
    except AttributeError:
        count = bin(num).count('1')

    # Part 2: Find a1 and b1 for the largest term in the sequence.
    # The largest term corresponds to the largest power of 2 in the binary expansion of N=10^100.
    # The exponent of this power of 2 is k_max = floor(log2(10^100)).
    # N.bit_length() - 1 is an efficient way to compute floor(log2(N)) for a positive integer N.
    N = 10**100
    k_max = N.bit_length() - 1

    # Now, we need to find a1 and b1 such that k_max = 2^(a1-1) + b1.
    # From our analysis, a1 - 1 = floor(log2(k_max)).
    # k_max.bit_length() - 1 is an efficient way to compute floor(log2(k_max)).
    if k_max == 0:
        a1_minus_1 = -1 # log2(0) is undefined, handle as special case if necessary
    else:
        a1_minus_1 = k_max.bit_length() - 1
    
    a1 = a1_minus_1 + 1
    
    # b1 = k_max - 2^(a1-1)
    b1 = k_max - (1 << a1_minus_1)

    # Print the final result in the specified format.
    print(f"{count} {a1} {b1}")

solve()

# The final answer will be printed by the function.
# The numeric result based on calculation is 115 9 76.
# Let's ensure the printout is just the numbers.

def final_solve():
    num_for_count = 5**100
    try:
        count = num_for_count.bit_count()
    except AttributeError:
        count = bin(num_for_count).count('1')
        
    N = 10**100
    k_max = N.bit_length() - 1
    
    a1_minus_1 = k_max.bit_length() - 1
    a1 = a1_minus_1 + 1
    b1 = k_max - (1 << a1_minus_1)
    
    # Final required output format
    print(count, a1, b1)

# Calling final_solve() would produce the desired output line.
# However, per instructions, I'll place the code directly in the block.
# final_solve() is essentially the same as solve() with less commentary.
