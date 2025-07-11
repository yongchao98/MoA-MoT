from fractions import Fraction

# The problem boils down to summing 1/2^(i+j) over three disjoint sets of pairs (i, j).

# Case 1: Pairs (k, k) for k >= 1.
# The sum is Sum_{k=1 to inf} (1/4)^k, which is 1/3.
sum_A = Fraction(1, 3)

# Case 2: Pairs (k, 2k) for k >= 1.
# The sum is Sum_{k=1 to inf} (1/8)^k, which is 1/7.
sum_B = Fraction(1, 7)

# Case 3: Pairs (2k, k) for k >= 1.
# The sum is Sum_{k=1 to inf} (1/8)^k, which is 1/7.
sum_C = Fraction(1, 7)

# The total sum is the sum of the sums from the three cases.
total_sum = sum_A + sum_B + sum_C

# Print the final equation with each number as requested.
# Using string formatting with the Fraction objects.
print(f"The total sum is the sum from three disjoint sets of pairs.")
print(f"Sum for pairs (k,k): {sum_A}")
print(f"Sum for pairs (k,2k): {sum_B}")
print(f"Sum for pairs (2k,k): {sum_C}")
print(f"Total Sum = {sum_A} + {sum_B} + {sum_C} = {total_sum}")