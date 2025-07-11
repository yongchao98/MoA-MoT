import math

# Given values
axial_tilt_epsilon = 23.5  # in degrees
period_precession = 26000 # in years
time_A_equator = -3000 # years from now
time_B_equator = 10000 # years from now

# Step 1: Calculate the time difference between the two events.
time_difference = time_B_equator - time_A_equator

# Step 2: Calculate the angle of precession corresponding to this time difference.
# The precession completes a full 360 degrees in T = 26000 years.
precession_angle = (time_difference / period_precession) * 360

# Step 3: Determine the angular separation of the North Celestial Pole (NCP) positions.
# The NCP precesses around the North Ecliptic Pole (NEP) in a circle of radius epsilon.
# The angular distance (d) between two points on this circle separated by a phase angle (theta)
# is given by the spherical law of cosines on the isosceles triangle NEP-NCP1-NCP2.
# cos(d) = cos(epsilon)^2 + sin(epsilon)^2 * cos(theta)
# When theta is 180 degrees (as calculated in Step 2), this simplifies to:
# cos(d) = cos(epsilon)^2 - sin(epsilon)^2 = cos(2 * epsilon)
# Therefore, the angular distance between the two NCP positions is d = 2 * epsilon.
angular_distance_NCPs = 2 * axial_tilt_epsilon

# Step 4: Relate the star positions to the NCP positions.
# A star is on the celestial equator when its angular distance from the NCP is 90 degrees.
# Star A is on the equator when the NCP is at position NCP_A.
# Star B is on the equator when the NCP is at position NCP_B.
# The problem's deep symmetry (swapping coordinates, and the T/2 time difference)
# implies that the fixed distance between Star A and Star B must be equal to the
# distance between the NCP positions NCP_A and NCP_B at the moments of their equator crossings.
angular_distance_stars = angular_distance_NCPs

# Print the final result and the equation
print("The angular distance between the two stars is derived from the separation of the North Celestial Pole positions at the times of the equator crossings.")
print(f"The time difference between crossings is {time_B_equator} - ({time_A_equator}) = {time_difference} years.")
print(f"This corresponds to a precession angle of ({time_difference}/{period_precession})*360 = {precession_angle} degrees.")
print("This means the NCP positions are diametrically opposite on the precession circle.")
print("The angular distance 'D' between these two NCP positions is given by D = 2 * epsilon.")
print(f"D = 2 * {axial_tilt_epsilon}")
print(f"D = {angular_distance_stars}")