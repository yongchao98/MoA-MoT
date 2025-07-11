import math

# Step 1: Define the values as fractions with small integers
density_num = 9
density_den = 10
radius_num = 1
radius_den = 2

# Step 2: Define constants from the volume formula
vol_const_num = 4
vol_const_den = 3

# Step 3: Use the simplified approximation for pi
pi_approx_num = 3
pi_approx_den = 1 # Representing 3 as 3/1

# The formula is: Mass = Density * Volume = Density * (4/3) * pi * radius^3
print("To estimate the mass, the parrot can use this equation:")
# The prompt requires printing each number in the equation.
print(f"Mass = ({density_num} / {density_den}) * ({vol_const_num} / {vol_const_den}) * ({pi_approx_num} / {pi_approx_den}) * ({radius_num} / {radius_den})^3 kg")

# Step 4: Perform the calculation
# First, calculate radius cubed
radius_cubed_num = radius_num ** 3
radius_cubed_den = radius_den ** 3

print("\nLet's calculate step by step:")
print(f"First, we calculate (1/2)^3 which is {radius_cubed_num}/{radius_cubed_den}.")
print(f"The equation becomes: Mass = {density_num}/{density_den} * {vol_const_num}/{vol_const_den} * {pi_approx_num}/{pi_approx_den} * {radius_cubed_num}/{radius_cubed_den} kg")

# Calculate the numerator and denominator of the final fraction
final_numerator = density_num * vol_const_num * pi_approx_num * radius_cubed_num
final_denominator = density_den * vol_const_den * pi_approx_den * radius_cubed_den

print("\nTo simplify, we can cancel out the two 3s:")
print(f"Mass = ({density_num}/{density_den}) * 4 * (1/{radius_cubed_den}) kg")
print(f"Mass = ({density_num}/{density_den}) * (4/{radius_cubed_den}) kg")
print(f"Mass = ({density_num}/{density_den}) * (1/2) kg")


# To get the final simplified fraction, find the greatest common divisor (GCD)
common_divisor = math.gcd(final_numerator, final_denominator)
simple_numerator = final_numerator // common_divisor
simple_denominator = final_denominator // common_divisor

print("\nThe final estimated mass is:")
print(f"Mass = {simple_numerator} / {simple_denominator} kg")
