import math

# This script finds a simple estimation for the mass of the rock
# and prints the calculation.

# Approximations for the formula: Mass = (4/3) * pi * r^3 * density
# where r = 1/2, pi ≈ 3, and density ≈ 1.

# The integers used for the calculation are all 4 or less.
vol_const_num = 4
vol_const_den = 3
pi_approx = 3
r_num = 1
r_den = 2
density_approx = 1

# Calculate the numerator and denominator of the result
numerator = vol_const_num * pi_approx * (r_num**3) * density_approx
denominator = vol_const_den * (r_den**3)

# Simplify the resulting fraction to its lowest terms
common_divisor = math.gcd(numerator, denominator)
final_num = numerator // common_divisor
final_den = denominator // common_divisor

# Print the final equation, showing each number used
print("The parrot can use this calculation:")
print(f"({vol_const_num} / {vol_const_den}) * {pi_approx} * ({r_num} / {r_den})^3 * {density_approx} = {final_num} / {final_den}")