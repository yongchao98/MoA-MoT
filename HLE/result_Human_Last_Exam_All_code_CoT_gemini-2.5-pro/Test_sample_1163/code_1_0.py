import math

# Earth's axial tilt in degrees, which is equal to the ecliptic latitude of the stars.
beta = 23.5

# The angular separation formula simplifies to 180 - 2*beta
two_beta = 2 * beta
angular_distance = 180 - two_beta

# Output the final calculation steps as requested.
# The equation for the angular distance is D = 180 - 2 * beta.
# So, D = 180 - 2 * 23.5 = 180 - 47 = 133.
print(f"The calculation for the angular distance is:")
print(f"180 - (2 * {beta}) = 180 - {two_beta} = {angular_distance}")
print(f"The angular distance between the two stars is {angular_distance} degrees.")