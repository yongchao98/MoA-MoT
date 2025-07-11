import math

# Step 1: Define the given constants.
# Earth's axial tilt in degrees. This is the value for epsilon (ε).
# From the problem, this is also the magnitude of the ecliptic latitude of the stars (|β|).
axial_tilt_deg = 23.5

# Step 2: The logic derived from the problem statement leads to a simple formula for the angular distance.
# The ecliptic latitudes are equal in magnitude to the axial tilt: |β| = ε = 23.5 degrees.
# The difference in ecliptic longitudes is 180 degrees.
# The angular distance γ is given by the formula γ = 180° - 2 * |β|.

# Step 3: Calculate the angular distance.
# First, calculate 2 * |β|
two_beta_deg = 2 * axial_tilt_deg

# Then, calculate the final angular distance γ in degrees.
angular_distance_deg = 180 - two_beta_deg

# Step 4: Print the final equation and the result.
# The formula is γ = 180 - 2 * |β|
print("The angular distance (γ) is calculated using the relationship derived from the stars' motion and positions.")
print("The magnitude of the stars' common ecliptic latitude |β| is equal to the Earth's axial tilt ε.")
print(f"|β| = ε = {axial_tilt_deg} degrees.")
print("The angular separation is given by the formula: γ = 180° - 2 * |β|")
print(f"γ = 180 - 2 * {axial_tilt_deg}")
print(f"γ = 180 - {two_beta_deg}")
print(f"γ = {angular_distance_deg} degrees")
