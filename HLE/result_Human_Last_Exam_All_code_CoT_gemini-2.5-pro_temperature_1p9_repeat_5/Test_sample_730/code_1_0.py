import math

# This script formulates the calculation for the parrot.
# It does not perform the calculation itself, but presents the equation
# with integers that meet the parrot's constraints.

# Given values as fractions
density = "9/10"
volume_constant = "4/3"
# The best approximation for pi that allows the calculation to succeed
pi_approx = "10/3" 
# The cube of the radius (1/2)^3
radius_cubed = "1/8"

# Construct the equation string
equation = f"Mass = {density} * {volume_constant} * {pi_approx} * {radius_cubed}"

# Final output
print("Yes, the parrot can estimate the mass with the required accuracy.")
print("The calculation it can perform is:")
print(equation)