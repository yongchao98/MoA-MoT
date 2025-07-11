import math

# The problem simplifies to the calculation of a definite integral. The analytical
# solution of the integral is given by the expression: 3 * ln(3 / (e^2 + e + 1)).
# This script computes the numerical value of this expression.

# Mathematical constant e
e = math.e

# Define the components of the expression
numerator = 3.0
e_squared = e**2
denominator = e_squared + e + 1.0

# Calculate the final result
result = 3 * math.log(numerator / denominator)

# As requested, output the final equation showing the numerical values.
print(f"The analytical result for the integral is: 3 * ln(3 / (e^2 + e + 1))")
print(f"Plugging in the numbers, the final calculation is:")
print(f"3 * ln({numerator} / ({e_squared} + {e} + 1.0)) = {result}")
