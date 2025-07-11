from fractions import Fraction

# The problem asks for the value of the sum: Sum = sum over (n, m) in S of 1/(4^n).
#
# Our analysis shows that for each integer n >= 1, there are exactly two
# possible values for m. The pairs (n, m) in S are (n, 0) and
# (n, (3^n + (-1)^n)/2).
#
# So, for each n, we add two terms of 1/(4^n) to the total sum.
# The sum can be written as:
# Sum = sum_{n=1 to infinity} 2 * (1/4)^n
#
# This is a geometric series. We can calculate the sum using the formula S = a / (1-r)
# for the series starting at n=1, where a is the first term and r is the common ratio.
# The sum of (1/4)^n for n from 1 to infinity is (1/4) / (1 - 1/4).
# We then multiply the result by 2.

# Let's define the components of the calculation.
factor = 2
ratio = Fraction(1, 4)

# Calculate the sum of the geometric series part: Sum_series = r + r^2 + r^3 + ...
series_sum_numerator = ratio
series_sum_denominator = 1 - ratio
series_sum = series_sum_numerator / series_sum_denominator

# The final answer is the factor multiplied by the series sum.
final_answer = factor * series_sum

print("The calculation for the sum is based on the geometric series formula.")
print("The sum is given by the expression: 2 * sum_{n=1 to infinity} (1/4)^n")
print(f"This can be written as: {factor} * ( (1/4) + (1/16) + (1/64) + ... )")
print("Let r = 1/4. The sum of the series part is given by the formula: r / (1 - r).")
print(f"The calculation is: {factor} * ( ({ratio.numerator}/{ratio.denominator}) / (1 - {ratio.numerator}/{ratio.denominator}) )")
print(f"= {factor} * ( ({ratio.numerator}/{ratio.denominator}) / ({series_sum_denominator.numerator}/{series_sum_denominator.denominator}) )")
print(f"= {factor} * ({series_sum.numerator}/{series_sum.denominator})")
print(f"= {final_answer.numerator}/{final_answer.denominator}")

final_value = float(final_answer)
print(f"\nThe final numeric value of the sum is {final_value}.")
