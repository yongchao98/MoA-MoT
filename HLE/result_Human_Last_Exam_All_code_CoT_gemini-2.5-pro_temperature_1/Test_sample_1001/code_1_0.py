from fractions import Fraction

# We identified the set S as the union of three disjoint families of pairs:
# 1. (k, k) for k >= 1
# 2. (k, 2k) for k >= 1
# 3. (2k, k) for k >= 1
# We calculate the sum of 1/2^(i+j) for each family and add them up.

# Sum for pairs (k, k) corresponds to the geometric series sum( (1/4)^k ) for k=1 to inf.
# Sum = (1/4) / (1 - 1/4) = 1/3.
sum1 = Fraction(1, 3)

# Sum for pairs (k, 2k) corresponds to the geometric series sum( (1/8)^k ) for k=1 to inf.
# Sum = (1/8) / (1 - 1/8) = 1/7.
sum2 = Fraction(1, 7)

# Sum for pairs (2k, k) is the same as for (k, 2k).
# Sum = 1/7.
sum3 = Fraction(1, 7)

# The total sum is the sum of the three parts.
total_sum = sum1 + sum2 + sum3

# Print the final equation with each number.
print(f"The total sum is calculated as the sum of three series:")
print(f"Sum for pairs (k, k): {sum1}")
print(f"Sum for pairs (k, 2k): {sum2}")
print(f"Sum for pairs (2k, k): {sum3}")
print(f"Final Equation: {sum1} + {sum2} + {sum3} = {total_sum}")
