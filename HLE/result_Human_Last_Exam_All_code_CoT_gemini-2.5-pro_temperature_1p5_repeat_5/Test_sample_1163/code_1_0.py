import math

# Given values
# Axial tilt in degrees
epsilon_deg = 23.5

# The period of precession in years
T = 26000
# Time when star A was on the celestial equator
t_A = -3000
# Time when star B will be on the celestial equator
t_B = 10000

# The problem's geometry leads to the conclusion that the angular distance 
# between the two stars is equal to twice the Earth's axial tilt.
# This arises from the symmetrical nature of the setup. The time difference 
# between the equator crossings is exactly half the precession period (13000 years),
# which places the North Celestial Pole at two diametrically opposite positions 
# (C_A and C_B) on its precession circle around the North Ecliptic Pole (NEP).
# The condition that the stars' coordinates will swap implies they have the same
# ecliptic latitude, meaning they are equidistant from the NEP.
# The conditions d(Star_A, C_A) = 90 degrees and d(Star_B, C_B) = 90 degrees,
# combined with the other symmetries, fix the angular distance between the stars
# to be equal to the angular distance between C_A and C_B, which is 2 * epsilon.

# Calculate the angular distance
angular_distance = 2 * epsilon_deg

# Output the equation and the result
print(f"The angular distance between the two stars can be calculated as:")
print(f"Distance = 2 * Axial Tilt")
print(f"Distance = 2 * {epsilon_deg}°")
print(f"Distance = {angular_distance}°")
