from fractions import Fraction

# The problem simplifies to calculating the sum of a geometric series.
# The set S contains two pairs for each n >= 1, so we sum 1/4^n twice for each n.
# This gives the formula: 2 * sum_{n=1 to infinity} (1/4)^n.

# The sum of an infinite geometric series sum_{k=1 to inf} r^k is r / (1-r).
# Here, r = 1/4.

coefficient = 2
ratio = Fraction(1, 4)

# Calculate the sum of the geometric series part.
series_sum = ratio / (1 - ratio)

# Calculate the final total sum.
total_sum = coefficient * series_sum

print(f"The value of the sum is determined by the formula: {coefficient} * (sum of (1/4)^n for n from 1 to infinity)")
print(f"The sum of the geometric series is r / (1 - r) where r = {ratio}.")
print(f"The final calculation is: {coefficient} * ( {ratio} / (1 - {ratio}) ) = {total_sum}")

print("\nThe final numerical value is:")
print(total_sum)