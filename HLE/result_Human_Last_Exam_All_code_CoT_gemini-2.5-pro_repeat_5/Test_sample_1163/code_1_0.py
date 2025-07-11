import math

# Given values
# Axial tilt of the Earth in degrees
axial_tilt_epsilon = 23.5

# The angular distance between the two stars is twice the axial tilt.
# This is because the stars are located at ecliptic latitudes of +epsilon and -epsilon
# and on the same ecliptic longitude.
angular_distance = 2 * axial_tilt_epsilon

# Print the equation with the numbers used
print(f"The final equation is: Angular Distance = 2 * {axial_tilt_epsilon}")

# Print the final result
print(f"The angular distance between the two stars is {angular_distance} degrees.")
