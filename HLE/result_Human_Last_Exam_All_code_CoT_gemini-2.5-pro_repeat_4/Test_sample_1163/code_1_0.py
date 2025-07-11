import math

# Given values
axial_tilt_deg = 23.5 # degrees (epsilon)
period_years = 26000 # years
time_A_years = -3000 # years ago
time_B_years = 10000 # years from now

# The reasoning leads to the conclusion that the angular distance between the stars
# is equal to the angular distance between the positions of the celestial pole
# at the times of the respective equator crossings.

# The time difference between the two events is 10000 - (-3000) = 13000 years.
time_diff_years = time_B_years - time_A_years

# The angular speed of precession in degrees per year.
omega_deg_per_year = 360.0 / period_years

# The angle subtended by the arc connecting the two pole positions at the center of the precession circle (the ecliptic pole).
pole_angle_diff_deg = omega_deg_per_year * time_diff_years

# The distance between the two pole positions P_A and P_B is the length of a chord on the celestial sphere
# connecting two points on a small circle of radius epsilon.
# The formula is 2 * epsilon * sin(pole_angle_diff_deg / 2).
# However, this is the straight-line distance in 3D space. The question asks for angular distance on the sphere.
# The angular distance is found using the law of cosines on the isosceles spherical triangle KP_A P_B, where K is the ecliptic pole.
# cos(D_poles) = cos(epsilon)*cos(epsilon) + sin(epsilon)*sin(epsilon)*cos(pole_angle_diff_deg)
# cos(D_poles) = cos^2(epsilon) + sin^2(epsilon)*cos(180)
# cos(D_poles) = cos^2(epsilon) - sin^2(epsilon) = cos(2*epsilon)
# So, the angular distance between the pole positions is D_poles = 2 * epsilon.

# Let's calculate this value.
angular_distance_deg = 2 * axial_tilt_deg

print(f"The axial tilt (epsilon) is {axial_tilt_deg} degrees.")
print(f"The angular distance between the pole positions at the time of the equator crossings is 2 * epsilon.")
print(f"Final Equation: Angular Distance = 2 * {axial_tilt_deg}")
print(f"Result: {angular_distance_deg}")
