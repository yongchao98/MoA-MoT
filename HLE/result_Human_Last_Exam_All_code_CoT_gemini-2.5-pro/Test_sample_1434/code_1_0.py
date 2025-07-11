from fractions import Fraction

# The problem is to calculate the sum S = sum_{(n, m) in S} 1/4^n.
#
# Our derivation shows that for each integer n >= 1, there are exactly two
# pairs (n, m) in the set S. These pairs are (n, 0) and (n, (3^n + (-1)^n)/2).
#
# Since for each n, there are two pairs in S, the term 1/4^n is contributed twice to the sum.
# The total sum is S = sum_{n=1 to infinity} (1/4^n + 1/4^n)
# This simplifies to S = sum_{n=1 to infinity} 2 * (1/4)^n
# We can factor out the constant 2: S = 2 * sum_{n=1 to infinity} (1/4)^n.
#
# The sum part is a geometric series sum_{n=1 to inf} r^n, which converges to r / (1 - r) for |r| < 1.
# In our case, the ratio r is 1/4.

# Define the coefficient and the ratio of the series
coefficient = 2
ratio = Fraction(1, 4)

# Calculate the sum of the geometric series: ratio / (1 - ratio)
# Since the sum starts from n=1, the first term 'a' is equal to the ratio 'r'.
series_sum_numerator = ratio
series_sum_denominator = 1 - ratio
series_sum = series_sum_numerator / series_sum_denominator

# Calculate the final sum by multiplying by the coefficient.
total_sum = coefficient * series_sum

# Output each number in the final equation as it is being calculated.
print(f"The total sum is given by the expression: {coefficient} * (sum from n=1 to infinity of ({ratio})^n)")
print(f"The sum of the geometric series is (a / (1 - r)) where a=r={ratio}:")
print(f"Numerator a = {series_sum_numerator}")
print(f"Denominator (1-r) = 1 - {ratio} = {series_sum_denominator}")
print(f"Sum of series = {series_sum_numerator} / {series_sum_denominator} = {series_sum}")
print(f"The final value is {coefficient} * {series_sum} = {total_sum}")