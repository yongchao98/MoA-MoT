import math

# Given parameters
axial_tilt_epsilon = 23.5  # degrees

# From the problem's constraints, we deduced that the ecliptic latitude (beta) of both stars
# must be equal in magnitude to the Earth's axial tilt.
# |beta| = epsilon
beta_deg = -axial_tilt_epsilon

# The angular separation (delta_sigma) between two points on a sphere with the same latitude (beta)
# and a longitude difference of 180 degrees is given by the formula:
# cos(delta_sigma) = -cos(2 * beta)
# This simplifies to delta_sigma = 180 - |2 * beta| in degrees.

# Step 1: Calculate 2 * beta
two_beta_deg = 2 * beta_deg

# Step 2: Calculate the absolute value of 2 * beta
abs_two_beta_deg = abs(two_beta_deg)

# Step 3: Calculate the angular distance
angular_distance_deg = 180 - abs_two_beta_deg

# Output the equation and the final result
# The problem asks for the equation to be printed.
print(f"The ecliptic latitude of the stars, |β|, is determined to be equal to the axial tilt, ε.")
print(f"|β| = {axial_tilt_epsilon}°")
print(f"The angular distance Δσ is calculated using the formula: Δσ = 180° - |2 * β|")
print(f"Δσ = 180° - |2 * ({beta_deg}°)|")
print(f"Δσ = 180° - |{two_beta_deg}°|")
print(f"Δσ = 180° - {abs_two_beta_deg}°")
print(f"Δσ = {angular_distance_deg}°")
