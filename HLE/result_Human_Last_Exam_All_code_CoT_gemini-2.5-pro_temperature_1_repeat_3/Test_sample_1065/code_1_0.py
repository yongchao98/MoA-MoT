import math

# The problem asks for the coordinates of the center of mass of a uniform
# quarter-circular arc of radius R, located in the first quadrant.
# The formula for the center of mass (x_cm, y_cm) of such an arc is:
# x_cm = 2*R/pi
# y_cm = 2*R/pi
#
# The user wants a "raw number", which implies a result independent of R.
# We provide the dimensionless coordinates x_cm/R and y_cm/R.
# This corresponds to calculating 2/pi.

# The numerator in the final equation for the coordinate ratio
numerator = 2

# The denominator in the final equation for the coordinate ratio
# is pi. We use the math library for an accurate value.
denominator = math.pi

# Calculate the value of the coordinate.
# Both horizontal (x) and vertical (y) coordinates are the same by symmetry.
coordinate_value = numerator / denominator

# Print the horizontal and vertical coordinates, separated by a comma.
# As requested, we output each number that makes up the final result.
# In this case, the horizontal and vertical coordinates are the numbers.
print(f"{coordinate_value},{coordinate_value}")