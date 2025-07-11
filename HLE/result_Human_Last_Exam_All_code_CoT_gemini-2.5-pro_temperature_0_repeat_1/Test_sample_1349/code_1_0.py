import math

# The supremum of X is given by the expression 40 / (36 * pi^2 + 25).
# Here we calculate its numerical value.

pi = math.pi

# Define the numerator and the denominator of the expression
numerator = 40
denominator_term1 = 36 * pi**2
denominator_term2 = 25
denominator = denominator_term1 + denominator_term2

# Calculate the final value
supremum_X = numerator / denominator

# Print the components of the final expression and the result
print("The supremum of X is given by the expression: numerator / (36 * pi^2 + 25)")
print(f"Numerator: {numerator}")
print(f"Denominator: 36 * {pi**2:.4f} + 25 = {denominator_term1:.4f} + {denominator_term2} = {denominator:.4f}")
print(f"The numerical value for the supremum of X is: {supremum_X}")

# Final answer format
# print(f'<<<{supremum_X}>>>')