from fractions import Fraction

# The problem reduces to calculating the sum of a geometric series.
# The sum is Sum = 2 * sum_{n=1 to infinity} (1/4)^n.
# The formula for an infinite geometric series sum_{n=1 to infinity} r^n is r / (1 - r).

# Define the numbers in the final equation
factor = 2
ratio_numerator = 1
ratio_denominator = 4

# The equation for the sum is:
# factor * ( (ratio_numerator / ratio_denominator) / (1 - ratio_numerator / ratio_denominator) )

# Calculate the sum using floating-point numbers
r = ratio_numerator / ratio_denominator
series_sum_float = r / (1 - r)
total_sum_float = factor * series_sum_float

print(f"The sum is given by the formula: {factor} * ( ({ratio_numerator}/{ratio_denominator}) / (1 - {ratio_numerator}/{ratio_denominator}) )")
print(f"The value of the sum is: {total_sum_float}")

# Calculate the sum using exact fractions for precision
r_frac = Fraction(ratio_numerator, ratio_denominator)
series_sum_frac = r_frac / (1 - r_frac)
total_sum_frac = factor * series_sum_frac

print(f"The exact value of the sum is: {total_sum_frac}")
