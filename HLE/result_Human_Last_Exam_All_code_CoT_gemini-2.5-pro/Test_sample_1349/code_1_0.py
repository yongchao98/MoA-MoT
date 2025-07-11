import math

# The supremum of X was found to be 40 / (25 + 36 * pi^2).
# This script calculates the numerical value of this expression.

# Define the numbers present in the final equation
numerator = 40
constant_term = 25
pi_coefficient = 36

# Calculate pi squared
pi_squared = math.pi**2

# Perform the calculation
denominator = constant_term + pi_coefficient * pi_squared
supremum_X = numerator / denominator

# Output each number in the final equation and the final result
print(f"The supremum of X is given by the expression: {numerator} / ({constant_term} + {pi_coefficient} * \u03c0^2)")
print(f"Substituting the value of \u03c0^2 \u2248 {pi_squared}, the expression becomes:")
print(f"{numerator} / ({constant_term} + {pi_coefficient} * {pi_squared})")
print(f"The numerical value of the supremum of X is: {supremum_X}")