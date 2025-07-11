import math

# Given values
axial_tilt_epsilon = 23.5  # in degrees

# From the derivation, the ecliptic latitude of the stars is the negative of the axial tilt.
beta_deg = -axial_tilt_epsilon

# The angular distance psi between two points with the same latitude beta
# and longitudes separated by 180 degrees is given by:
# cos(psi) = -cos(2*beta)
# Therefore, psi = 180 - 2*|beta|
# Since beta = -epsilon, |beta| = epsilon.
angular_distance_deg = 180 - 2 * axial_tilt_epsilon

# The equation is: angular distance = 180 - 2 * 23.5
part1 = 180
part2 = 2
part3 = 23.5
result = part1 - part2 * part3

print(f"The angular distance between the two stars is derived from the equation:")
print(f"{part1} - {part2} * {part3} = {result}")
print(f"The angular distance between the two stars is {result} degrees.")
