# The Earth's axial tilt (obliquity of the ecliptic) in degrees.
axial_tilt_epsilon = 23.5

# The problem states that star A was on the celestial equator at a past time (t_A)
# and star B will be on the celestial equator at a future time (t_B).
# The time difference between these events is 13000 years, which is half of the
# Earth's precession period of 26000 years.
# This means the North Celestial Pole's (NCP) position at time t_A, let's call it P_A,
# is diametrically opposite to its position at time t_B, let's call it P_B,
# with respect to the North Ecliptic Pole (NEP).

# The angular distance between P_A and P_B on the celestial sphere is the
# diameter of the precession circle.
# Distance(P_A, P_B) = 2 * radius_of_precession_circle
# The radius of the precession circle is the axial tilt, epsilon.
distance_between_poles = 2 * axial_tilt_epsilon

# The problem states that the stars will swap their equatorial coordinates at a future time.
# This implies a strong symmetry in their positions.
# Let S_A and S_B be the positions of the stars.
# The condition that S_A is on the equator at time t_A means its angular distance from P_A is 90 degrees.
# The condition that S_B is on the equator at time t_B means its angular distance from P_B is 90 degrees.
# Due to the inherent symmetry of the problem, the geometric relationship between
# star S_A and pole P_A is the same as the relationship between star S_B and pole P_B.
# This implies that the angular distance between the stars S_A and S_B is equal
# to the angular distance between the pole positions P_A and P_B.

angular_distance = distance_between_poles

print("The logic of the solution is as follows:")
print("1. The time between Star A being on the equator and Star B being on the equator is 13000 years, which is half the precession period (26000 years).")
print("2. This means the positions of the North Celestial Pole (NCP) at these two moments are diametrically opposite on the precession circle.")
print("3. The radius of the precession circle is the axial tilt, epsilon = 23.5 degrees.")
print("4. The angular distance between these two opposite NCP positions is the diameter of this circle, which is 2 * epsilon.")
print("5. Due to the coordinate-swapping symmetry, the angular distance between the two stars must equal the angular distance between their corresponding NCP positions.")
print("6. Therefore, the angular distance between the stars is 2 * epsilon.")
print("\nThe calculation is:")
print("Angular Distance = 2 * Axial Tilt")
print("Angular Distance = 2 * " + str(axial_tilt_epsilon))
print("Angular Distance = " + str(angular_distance) + " degrees")
print("\nThe final answer for the angular distance is " + str(angular_distance) + " degrees.")