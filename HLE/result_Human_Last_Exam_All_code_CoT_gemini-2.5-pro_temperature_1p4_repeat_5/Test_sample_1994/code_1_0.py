import math

# The number of sequences is 1, as the representation of 10^100 as a sum
# of distinct powers of 2 (each mappable to a unique tetration term) is unique.
count = 1

# Find the exponent k for the largest power of 2 less than or equal to 10^100
# k = floor(log2(10^100)) = floor(100 * log2(10))
log2_of_N = 100 * math.log2(10)
k1 = math.floor(log2_of_N)

# The largest term in the sequence is 2^k1. We need to find the values
# a1 and b1 such that tet(2, a1) * pow(2, b1) = 2^k1.
# This means we must solve for a1 and b1 in the equation:
# k1 = 2^(a1 - 1) + b1, with the constraint b1 < 2^(a1 - 1).

# To find a1, we find the largest power of 2 less than or equal to k1.
# The exponent of that power of 2 will be (a1 - 1).
# a1 - 1 = floor(log2(k1))
a1_minus_1 = math.floor(math.log2(k1))
a1 = a1_minus_1 + 1

# Now we can calculate b1.
# b1 = k1 - 2^(a1 - 1)
b1 = k1 - (2**a1_minus_1)

# Print the results in the required format.
print(f"{count} {a1} {b1}")
