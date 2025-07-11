from fractions import Fraction

# The problem's solution reduces to calculating the sum of a geometric series.
# The sum is given by the formula: 2 * sum_{n=1 to infinity} (1/4)^n.

# Define the constants for the calculation.
factor = 2
ratio = Fraction(1, 4)

# Calculate the sum of the infinite geometric series sum_{n=1 to inf} r^n = r / (1-r).
geometric_series_sum = ratio / (1 - ratio)

# Calculate the final total sum.
total_sum = factor * geometric_series_sum

# Output the equation with each number substituted, showing the steps of calculation.
print(f"The formula for the sum is: {factor} * ( Sum_{{n=1}}^{{\infty}} ( {ratio} )^n )")
print(f"This is a geometric series with sum equal to: {factor} * ( ({ratio}) / (1 - {ratio}) )")
print(f"Calculating the denominator: 1 - {ratio} = {1-ratio}")
print(f"Calculating the series sum: ({ratio}) / ({1-ratio}) = {geometric_series_sum}")
print(f"The final calculation is: {factor} * {geometric_series_sum} = {total_sum}")
print(f"\nThe value of the sum is {total_sum}.")
