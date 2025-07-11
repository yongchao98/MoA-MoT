import math

# The problem asks for the limit of n * P(n) as n -> infinity.
# Our derivation shows that this limit is equal to the constant value (2 * sqrt(3)) / pi.

# Value of the numerator coefficient
numerator_coeff = 2
# Value inside the square root in the numerator
numerator_sqrt_val = 3

# Value of the denominator
denominator_val = math.pi

# We will calculate the final numerical value and print the equation.
numerator = numerator_coeff * math.sqrt(numerator_sqrt_val)
result = numerator / denominator_val

print("The limit of n*P(n) is given by the expression: (2 * sqrt(3)) / pi")
print(f"Let's break down the calculation:")
print(f"Numerator: {numerator_coeff} * sqrt({numerator_sqrt_val}) = {numerator}")
print(f"Denominator: pi = {denominator_val}")
print(f"Final Result: {numerator} / {denominator_val} = {result}")
