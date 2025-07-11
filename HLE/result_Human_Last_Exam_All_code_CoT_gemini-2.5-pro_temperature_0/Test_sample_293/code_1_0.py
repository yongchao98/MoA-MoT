import math

# The integral has been symbolically evaluated to -pi/2.
# This script calculates its numerical value.

# Define the components of the final expression
sign_char = "-"
numerator_symbol = "pi"
denominator_value = 2
numerator_value = math.pi

# Calculate the final result
result = -numerator_value / denominator_value

# Output the final equation and its numerical value,
# including each component as requested.
print(f"The integral evaluates to the expression: {sign_char}{numerator_symbol} / {denominator_value}")
print(f"The components are:")
print(f"  - Sign: {sign_char}")
print(f"  - Numerator ({numerator_symbol}): {numerator_value}")
print(f"  - Denominator: {denominator_value}")
print(f"\nThe final value of the integral is: {result}")