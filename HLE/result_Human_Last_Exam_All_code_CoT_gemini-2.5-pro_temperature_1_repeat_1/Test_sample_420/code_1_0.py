import math

# Based on the geometric and calculus-based derivation, the final ratio is an exact expression.
# The problem asks to output the numbers in the final equation.
# The exact formula for the ratio is (9 + pi - 3*sqrt(3)) / 12.

# Define the numbers in the final equation's numerator and denominator.
numerator_constant_term = 9
numerator_pi_coefficient = 1
numerator_sqrt_coefficient = -3
value_inside_sqrt = 3
denominator = 12

# Print the final equation in its exact form.
print("The exact formula for the ratio of the area of D to the area of S is:")
print(f"({numerator_constant_term} + {numerator_pi_coefficient}*pi + ({numerator_sqrt_coefficient})*sqrt({value_inside_sqrt})) / {denominator}")

# For reference, we can also print the numerical approximation.
numerical_value = (9 + math.pi - 3 * math.sqrt(3)) / 12
print("\nThe approximate numerical value is:")
print(numerical_value)