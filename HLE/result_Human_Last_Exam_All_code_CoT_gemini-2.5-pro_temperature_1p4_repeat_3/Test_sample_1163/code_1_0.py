import math

# Given parameters
# Axial tilt of the Earth in degrees
epsilon_deg = 23.5

# The derivation shows that the two stars have the same ecliptic latitude,
# equal to the Earth's axial tilt (beta = epsilon), and their ecliptic
# longitudes are 180 degrees apart.

# The formula for the angular distance 'theta' between them is:
# theta = 180 - 2 * epsilon

# Step 1: Calculate 2 * epsilon
two_epsilon_deg = 2 * epsilon_deg

# Step 2: Calculate the final angular distance
angular_distance_deg = 180 - two_epsilon_deg

# Output the steps of the final calculation
print("The angular distance, theta, is derived from the formula: theta = 180 - 2 * epsilon")
print(f"Given epsilon = {epsilon_deg} degrees.")
print(f"First, calculate 2 * epsilon: 2 * {epsilon_deg} = {two_epsilon_deg} degrees.")
print(f"Then, calculate theta: 180 - {two_epsilon_deg} = {angular_distance_deg} degrees.")
print("\nThe final angular distance between the two stars is:")
print(f"{angular_distance_deg} degrees")
