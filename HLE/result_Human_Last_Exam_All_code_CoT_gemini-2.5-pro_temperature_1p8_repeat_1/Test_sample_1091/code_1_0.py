import math

# The problem is to find the limit of n * P(n) as n -> infinity.
# Through analytical derivation, the limit is found to be the expression 2*sqrt(3)/pi.

# Define the numbers in the final expression
numerator_factor = 2
sqrt_argument = 3
denominator_value = math.pi

# Calculate the result
limit_value = numerator_factor * math.sqrt(sqrt_argument) / denominator_value

# Output the final equation and the result.
print(f"The limit is given by the expression: {numerator_factor} * sqrt({sqrt_argument}) / pi")
print(f"Numerical value of the limit: {limit_value}")
