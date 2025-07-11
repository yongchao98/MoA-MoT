import math

# This script generates the calculation steps for the parrot.

# --- Problem values as fractions ---
# Density = 0.9 = 9/10
density_num = 9
density_den = 10

# Radius = 0.5 = 1/2
radius_num = 1
radius_den = 2

# --- Constants from the volume formula (4/3) * pi * r^3 ---
vol_const_num = 4
vol_const_den = 3

# --- Parrot-friendly approximation for pi ---
# Based on analysis, pi = 3 is a good approximation.
# It results in an error of ~4.5%, which is less than the 10% limit.
# The integers (3 and 1) are small.
pi_approx_num = 3
pi_approx_den = 1

# --- Construct and print the final equation for the parrot ---
# The parrot must see all the numbers involved in the calculation.
print("Hello Parrot! To find the mass of the rock, please calculate the following:")
print()
print(f"Mass = ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_approx_num}/{pi_approx_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den}) * ({radius_num}/{radius_den})")
print()

# --- Calculate and print the simplified final answer ---
# Calculate the numerator and denominator of the final fraction
final_num = density_num * vol_const_num * pi_approx_num * (radius_num**3)
final_den = density_den * vol_const_den * pi_approx_den * (radius_den**3)

# Simplify the fraction by dividing by the greatest common divisor (GCD)
common_divisor = math.gcd(final_num, final_den)
simple_num = final_num // common_divisor
simple_den = final_den // common_divisor

print(f"The final simplified answer is: {simple_num}/{simple_den} kg")
