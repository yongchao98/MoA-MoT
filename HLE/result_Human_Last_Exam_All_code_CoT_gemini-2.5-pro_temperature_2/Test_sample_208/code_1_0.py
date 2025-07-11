import math

# This script calculates the radius of a large circle required to tightly pack
# fourteen smaller circles of radius 1.

# The solution is based on the known optimal packing configuration for n=14.
# The formula for the radius 'R' of the containing circle is:
# R = 1 + 2 * sqrt((4 + 2*sqrt(3)) / 3)

print("The formula for the radius 'R' is: R = 1 + 2 * sqrt((4 + 2*sqrt(3)) / 3)")
print("This formula is used to find the radius of the large circle for packing 14 unit circles.")
print("\nLet's break down the calculation, showing each number in the equation:")

# Define the numbers from the formula
c1 = 1
c2_mult = 2
c4 = 4
c2_sqrt = 2
c3_inner = 3
c3_div = 3

# Display the numbers from the formula
print(f"The constant added at the end is: {c1}")
print(f"The multiplier outside the square root is: {c2_mult}")
print(f"The number inside the square root's numerator is: {c4} + {c2_sqrt}*sqrt({c3_inner})")
print(f"The denominator inside the square root is: {c3_div}")

# Perform the calculation
sqrt_3 = math.sqrt(c3_inner)
numerator = c4 + c2_sqrt * sqrt_3
fraction = numerator / c3_div
sqrt_of_fraction = math.sqrt(fraction)
final_radius = c1 + c2_mult * sqrt_of_fraction

# Print the final result
print(f"\nThe full calculated radius R is: {final_radius}")
print("\nThe problem asks for the radius up to 4 significant digits.")
print(f"The value rounded to 4 significant digits is: {final_radius:.4g}")
