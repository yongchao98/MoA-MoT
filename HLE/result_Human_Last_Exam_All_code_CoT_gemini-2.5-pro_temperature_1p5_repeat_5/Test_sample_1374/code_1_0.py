import math

# This script calculates the furthest distance from point A for a shape
# that maximizes the gravitational field at A, given a fixed volume of 1 m^3.

# The formula for this distance (R_max) is: R_max = (15 / (4 * pi))^(1/3)

# Here are the numbers used in that final equation:
numerator = 15
denominator_factor = 4
exponent_numerator = 1
exponent_denominator = 3

# Calculate the result
pi_value = math.pi
exponent = exponent_numerator / exponent_denominator
result = (numerator / (denominator_factor * pi_value)) ** exponent

# Output the equation and the final answer, as requested.
print(f"The calculation is based on the equation: ({numerator} / ({denominator_factor} * pi)) ^ ({exponent_numerator}/{exponent_denominator})")
print("\nThe furthest distance from point A is (in meters):")
print(result)