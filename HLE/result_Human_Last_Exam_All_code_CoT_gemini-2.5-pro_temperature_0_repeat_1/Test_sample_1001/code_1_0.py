from fractions import Fraction

# The problem reduces to finding three types of pairs (i, j) that satisfy the condition.
# The set S is the disjoint union of three sets of pairs:
# 1. (i, i) for i >= 1
# 2. (i, 2i) for i >= 1
# 3. (2i, i) for i >= 1

# We calculate the sum of 1/2^(i+j) over each set.

# For pairs (i, i), the sum is Sum_{i=1 to inf} 1/2^(2i) = Sum (1/4)^i
# This is a geometric series with a=1/4, r=1/4. Sum = a/(1-r) = (1/4)/(3/4) = 1/3.
sum_diag = Fraction(1, 3)

# For pairs (i, 2i), the sum is Sum_{i=1 to inf} 1/2^(3i) = Sum (1/8)^i
# This is a geometric series with a=1/8, r=1/8. Sum = a/(1-r) = (1/8)/(7/8) = 1/7.
sum_off_diag1 = Fraction(1, 7)

# For pairs (2i, i), the sum is Sum_{i=1 to inf} 1/2^(3i) = Sum (1/8)^i
# This is the same series, with sum 1/7.
sum_off_diag2 = Fraction(1, 7)

# The total sum is the sum of the sums over these three disjoint sets.
total_sum = sum_diag + sum_off_diag1 + sum_off_diag2

# Output the final calculation step-by-step
print("The total sum is the sum of three geometric series.")
print(f"Sum over pairs (i, i): {sum_diag.numerator}/{sum_diag.denominator}")
print(f"Sum over pairs (i, 2i): {sum_off_diag1.numerator}/{sum_off_diag1.denominator}")
print(f"Sum over pairs (2i, i): {sum_off_diag2.numerator}/{sum_off_diag2.denominator}")
print("\nThe final equation is:")
print(f"{sum_diag.numerator}/{sum_diag.denominator} + {sum_off_diag1.numerator}/{sum_off_diag1.denominator} + {sum_off_diag2.numerator}/{sum_off_diag2.denominator} = {total_sum.numerator}/{total_sum.denominator}")
