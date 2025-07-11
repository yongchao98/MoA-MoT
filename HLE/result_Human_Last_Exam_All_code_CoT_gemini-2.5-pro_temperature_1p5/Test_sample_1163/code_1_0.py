# Constants given in the problem
axial_tilt_epsilon = 23.5  # in degrees

# Step 1: Determine the distance of each star from the North Ecliptic Pole (NEP).
# Based on the symmetry argument, the stars and the relevant positions of the
# North Celestial Pole (P_A and P_B) lie on a single great circle.
# dist(Star, NCP) = 90 degrees for a star on the equator.
# The distance from the NCP position (e.g., P_A) to the NEP is the axial tilt (epsilon).
# The relationship on the great circle is: dist(Star, NEP) + dist(NEP, P_A) = dist(Star, P_A)
# So, dist_star_nep = 90 - epsilon.
dist_star_A_from_NEP = 90.0 - axial_tilt_epsilon
dist_star_B_from_NEP = 90.0 - axial_tilt_epsilon

# Step 2: Calculate the total angular distance between the two stars.
# The stars lie on the same great circle but on opposite sides of the NEP.
# So, the total distance is the sum of their individual distances from the NEP.
angular_distance = dist_star_A_from_NEP + dist_star_B_from_NEP

# Alternative calculation based on the derived formula
# angular_distance = (90 - epsilon) + (90 - epsilon) = 180 - 2 * epsilon
calculated_distance = 180.0 - 2 * axial_tilt_epsilon

# Print the final calculation step by step
print("The axial tilt (epsilon) is {} degrees.".format(axial_tilt_epsilon))
print("The distance of each star from the North Ecliptic Pole is 90 - epsilon.")
print("dist(Star A, NEP) = 90 - {} = {} degrees.".format(axial_tilt_epsilon, dist_star_A_from_NEP))
print("dist(Star B, NEP) = 90 - {} = {} degrees.".format(axial_tilt_epsilon, dist_star_B_from_NEP))
print("The total angular distance between the stars is the sum of these distances.")
print("Total Distance = {} + {} = {} degrees.".format(dist_star_A_from_NEP, dist_star_B_from_NEP, angular_distance))
print("Final Equation: 180 - 2 * {} = {}".format(axial_tilt_epsilon, calculated_distance))

# Final answer
print("\nThe angular distance between the two stars is {} degrees.".format(angular_distance))