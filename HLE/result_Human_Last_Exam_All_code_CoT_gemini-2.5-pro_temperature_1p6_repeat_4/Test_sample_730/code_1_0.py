import math

# Plan: The parrot needs to estimate the mass with an error of at most 10%.
# The mass is calculated by: mass = density * (4/3) * pi * radius^3.
# The parrot must use fractions with integers as small as possible.
# After searching for the optimal approximations for density (0.9) and pi (~3.14),
# the smallest possible 'largest integer' (z) in the calculation is 7.

# The calculation uses the following numbers:
# - Density (0.9) is approximated as 6/7.
# - The volume constant is 4/3.
# - Pi is approximated as 3 (or 3/1).
# - Radius (0.5) is 1/2.
# All integers used {6, 7, 4, 3, 3, 1, 2} are 7 or less.

# Numbers used by the parrot for the optimal calculation
density_num = 6
density_den = 7

sphere_const_num = 4
sphere_const_den = 3

pi_approx = 3

radius_num = 1
radius_den = 2

# Print the equation the parrot would use
print("Yes, the parrot can perform the estimation.")
print("The calculation that uses the smallest possible integers is:")
print(f"mass = ({density_num}/{density_den}) * ({sphere_const_num}/{sphere_const_den}) * {pi_approx} * ({radius_num}/{radius_den})^3")

# Calculate the final fraction from the parrot's equation
final_numerator = density_num * sphere_const_num * pi_approx * (radius_num**3)
final_denominator = density_den * sphere_const_den * (radius_den**3)

# Simplify the fraction for the final result
common_divisor = math.gcd(final_numerator, final_denominator)
simplified_num = final_numerator // common_divisor
simplified_den = final_denominator // common_divisor

# Print the final result of the equation
print(f"     = {simplified_num}/{simplified_den} kg")