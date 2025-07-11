import math

# The problem asks for the limit of n*P(n) as n goes to infinity.
# Based on the analysis using the Local Central Limit Theorem, this limit evaluates to the expression (2 * sqrt(3)) / pi.
# This script calculates and prints the value of this expression.

# The numbers that appear in the final derived mathematical expression for the limit.
numerator_factor = 2
sqrt_argument = 3
denominator_value = math.pi

# Calculate the final result.
result = (numerator_factor * math.sqrt(sqrt_argument)) / denominator_value

# Output the equation representing the limit and its numerical value.
# The following line prints each number used in the final equation.
print(f"Final equation for the limit: {numerator_factor} * sqrt({sqrt_argument}) / {denominator_value}")
print(f"Result: {result}")