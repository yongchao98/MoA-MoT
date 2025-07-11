import math
from fractions import Fraction

# This script calculates the exact ratio of the area of region D
# to the total surface area of the cube S.

# The final form of the ratio is A + B*pi + C*sqrt(3).
# Based on the derivation, the coefficients are:
# A = 3/4, B = 1/12, C = -1/4.

# We define the components of the final ratio.
rational_part = Fraction(3, 4)
pi_coefficient = Fraction(1, 12)
sqrt3_coefficient = Fraction(-1, 4)

# Extract numerators and denominators for printing the equation.
num_a, den_a = rational_part.numerator, rational_part.denominator
num_b, den_b = pi_coefficient.numerator, pi_coefficient.denominator
num_c, den_c = sqrt3_coefficient.numerator, sqrt3_coefficient.denominator

# Get the constant number for the square root.
sqrt_val = 3

print("The exact ratio Area(D) / Area(S) is calculated as:")
print(f"({num_a}/{den_a}) + ({num_b}/{den_b})*pi + ({num_c}/{den_c})*sqrt({sqrt_val})")
print("\nSo the final equation for the ratio is:")
print(f"{num_a}/{den_a} + pi/{den_b} - sqrt({sqrt_val})/{den_c}")

# Calculate the numerical value for reference.
numerical_value = float(rational_part) + float(pi_coefficient) * math.pi + float(sqrt3_coefficient) * math.sqrt(3)
print(f"\nFor reference, the approximate numerical value is: {numerical_value}")

# Return the exact final answer in the requested format.
final_answer_string = f"{num_a}/{den_a} + pi/{den_b} - sqrt({sqrt_val})/{abs(den_c)}"