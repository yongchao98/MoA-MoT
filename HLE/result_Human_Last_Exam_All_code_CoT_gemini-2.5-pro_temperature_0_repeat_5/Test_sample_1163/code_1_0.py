import math

# Given values
epsilon_deg = 23.5  # Earth's axial tilt in degrees

# The derivation shows that the angular distance 'd' can be found using the formula:
# d = 180 - 2 * epsilon
# This is because the ecliptic latitudes of the stars are both equal to epsilon,
# and their ecliptic longitudes differ by 180 degrees.

# Calculate 2 * epsilon
two_epsilon = 2 * epsilon_deg

# Calculate the final angular distance
angular_distance = 180 - two_epsilon

# Print the final equation with all the numbers
print("The final equation for the angular distance d is:")
print(f"d = 180 - 2 * {epsilon_deg}")
print(f"d = 180 - {two_epsilon}")
print(f"d = {angular_distance}")