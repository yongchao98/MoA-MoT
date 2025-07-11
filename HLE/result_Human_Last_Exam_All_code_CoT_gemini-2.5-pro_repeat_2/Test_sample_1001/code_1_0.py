from fractions import Fraction

# The problem boils down to summing three geometric series corresponding
# to three disjoint sets of pairs (i, j).

# 1. For pairs (k, k), the sum of 1/2^(i+j) is sum_{k=1 to inf} (1/4)^k
# which equals 1/3.
sum_for_S1 = Fraction(1, 3)

# 2. For pairs (k, 2k), the sum of 1/2^(i+j) is sum_{k=1 to inf} (1/8)^k
# which equals 1/7.
sum_for_S2 = Fraction(1, 7)

# 3. For pairs (2k, k), the sum is the same as for (k, 2k), which is 1/7.
sum_for_S3 = Fraction(1, 7)

# The total sum is the sum of these three parts.
total_sum = sum_for_S1 + sum_for_S2 + sum_for_S3

# Print the final equation with each component of the sum.
print(f"{sum_for_S1} + {sum_for_S2} + {sum_for_S3} = {total_sum}")