# Earth's axial tilt in degrees, given in the problem.
axial_tilt_epsilon = 23.5

# From the reasoning explained in the plan, the magnitude of the ecliptic latitude |beta|
# of both stars must be equal to the axial tilt.
beta = axial_tilt_epsilon

# The angular distance 'd' between the two stars, whose ecliptic longitudes differ by 180 degrees
# and have the same ecliptic latitude beta, is given by the formula: d = 180 - 2 * |beta|.
term_1 = 180
term_2 = 2
angular_distance = term_1 - term_2 * beta

# Print the final equation with all its components, as requested.
print(f"The calculation for the angular distance is:")
print(f"{term_1} - {term_2} * {beta} = {angular_distance}")
print(f"\nThe angular distance between the two stars is {angular_distance} degrees.")