import math

# The value of the integral I is given by the expression pi * ln(1 + sqrt(2)).
# This script calculates its numerical value.

# The numbers that make up the final equation are pi, 1, and 2.
pi = math.pi
one = 1.0
two = 2.0

# Perform the calculation.
sqrt_two = math.sqrt(two)
result = pi * math.log(one + sqrt_two)

# Output the equation with each number's value, as requested.
# The symbolic form is I = pi * ln(1 + sqrt(2)).
print(f"The integral is calculated using the formula: I = \u03c0 * ln(1 + \u221a2)")
print(f"I = {pi} * ln({one} + sqrt({two}))")

# Print the final numerical result.
print("The numerical value is:")
print(result)
