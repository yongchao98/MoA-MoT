# The axial tilt of the Earth in degrees.
axial_tilt_epsilon = 23.5

# Based on the analysis of the problem's constraints, the ecliptic latitude (beta)
# of both stars must be equal to the Earth's axial tilt.
beta = axial_tilt_epsilon

# The angular separation (delta_sigma) between two stars with the same latitude beta
# and longitudes separated by 180 degrees is given by the formula:
# delta_sigma = 180 - 2 * beta
# Here we calculate the values.
val_2_beta = 2 * beta
angular_distance = 180 - val_2_beta

# Output the final equation and the result
print(f"The angular distance between the two stars is calculated as:")
print(f"180 degrees - 2 * {beta} degrees = 180 degrees - {val_2_beta} degrees = {angular_distance} degrees")